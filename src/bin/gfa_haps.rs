// gfa_haps: extract two haplotype walks per connected component of a GFA.
//
// Strategy (two-stage):
//   1. For each connected component, generate K diverse candidate walks
//      via randomized heuristic search. Each walk traverses the bidirected
//      graph respecting orientation, prefers unvisited successors, and
//      minimizes within-walk revisits.
//   2. Select the best pair of walks (w1, w2) by a lexicographic objective:
//        a. maximize union coverage (bp covered by either walk)
//        b. maximize symmetric difference (bp covered by exactly one walk)
//           -- this drives the "two haplotypes diverge through bubbles" behavior
//        c. minimize total revisits across both walks
//
// Endpoints policy:
//   --endpoints tip|any
//     tip: walks must start and end at tips (nodes with no edges on the
//          relevant side). If a component has no tips, falls back to 'any'.
//     any: walks can start/end anywhere.
//
// Outputs (per input GFA file):
//   - <basename>.with_paths.gfa   : input + two P-lines per component
//   - <basename>.haps.fa          : one FASTA record per (component, hap)

use gfaphaser::{
    connected_components, entry_exit_for_orient, halfedge, parse_gfa, Edge, Graph, ParseOptions,
    Side, Step, Walk,
};
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process;

fn die(msg: &str) -> ! {
    eprintln!("error: {}", msg);
    process::exit(1);
}

// ------------ Simple deterministic PRNG (xorshift64*) ------------
//
// Avoids pulling in the `rand` crate so we stay dependency-free.

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        // Avoid 0; xorshift gets stuck.
        Rng(if seed == 0 { 0x9E3779B97F4A7C15 } else { seed })
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn gen_range(&mut self, lo: usize, hi: usize) -> usize {
        // hi exclusive
        let span = (hi - lo) as u64;
        lo + (self.next_u64() % span) as usize
    }
    /// Weighted pick. Weights must be > 0.
    fn weighted_pick(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        if total <= 0.0 {
            return self.gen_range(0, weights.len());
        }
        let mut r = (self.next_u64() as f64 / u64::MAX as f64) * total;
        for (i, &w) in weights.iter().enumerate() {
            r -= w;
            if r <= 0.0 {
                return i;
            }
        }
        weights.len() - 1
    }
}

// ------------ Candidate walk generation ------------

#[derive(Clone, Copy)]
enum EndpointMode {
    Tip,
    Any,
}

/// Find all (node, exit_side) candidates that look like good walk starts.
/// In "tip" mode we want nodes that have no edges on one side -- we'll
/// enter from the empty side (conceptually) and exit from the other side.
/// In "any" mode every (node, side) is a candidate.
fn walk_starts(g: &Graph, members: &[usize], mode: EndpointMode) -> Vec<(usize, Side)> {
    let mut starts = Vec::new();
    match mode {
        EndpointMode::Tip => {
            for &u in members {
                for side in [0u8, 1u8] {
                    // If side has no edges, we treat that as the "entry" side
                    // (the dead end) and we will exit out the *other* side.
                    if g.is_tip_on(u, side) && !g.is_tip_on(u, 1 - side) {
                        // exit through opposite side
                        starts.push((u, 1 - side));
                    }
                }
                // An isolated node (both sides empty) is also a valid start in tip mode
                if g.is_tip_on(u, 0) && g.is_tip_on(u, 1) {
                    starts.push((u, 1));
                }
            }
        }
        EndpointMode::Any => {
            for &u in members {
                starts.push((u, 0));
                starts.push((u, 1));
            }
        }
    }
    starts
}

/// Side encodes the *exit* side for the first node of the walk.
fn side_to_orient(exit_side: Side) -> char {
    // Exit side 1 = right => orientation '+'
    // Exit side 0 = left  => orientation '-'
    if exit_side == 1 {
        '+'
    } else {
        '-'
    }
}

/// Generate one randomized walk starting at (start_node, exit_side).
///
/// Termination rule (important): the walk stops when no edges are available
/// (dead end) OR when *all* available edges lead to already-visited nodes
/// AND the `allow_revisits` budget is exhausted. This prevents pathological
/// looping in cyclic components (e.g. bubbles).
///
/// `allow_revisits`: number of revisit steps the walk may take. Setting this
/// to a small value yields simple-path-like walks; larger values let walks
/// dip back into the graph if it can unlock new coverage downstream.
/// `global_avoid`: optional set of nodes used by another walk; the walk
/// strongly prefers neighbors NOT in this set (haplotype-divergence guide).
fn generate_walk(
    g: &Graph,
    start_node: usize,
    start_exit_side: Side,
    rng: &mut Rng,
    allow_revisits: usize,
    global_avoid: Option<&HashSet<usize>>,
    global_avoid_weight: f64,
) -> Walk {
    let mut walk = Walk::new();
    let mut visited_in_walk: HashSet<usize> = HashSet::new();

    let first_orient = side_to_orient(start_exit_side);
    walk.steps.push(Step {
        node: start_node,
        orient: first_orient,
        overlap_in: 0,
        cigar_in: String::new(),
    });
    visited_in_walk.insert(start_node);

    let mut cur = start_node;
    let mut cur_exit = start_exit_side;
    let mut revisits_used: usize = 0;
    // Hard ceiling to prevent any pathological case from running forever.
    let hard_cap = g.n().saturating_mul(8).max(64);

    while walk.steps.len() < hard_cap {
        let edges = g.neighbors(cur, cur_exit);
        if edges.is_empty() {
            break;
        }

        // Separate unvisited vs visited destinations.
        let unvisited: Vec<&Edge> = edges
            .iter()
            .filter(|e| !visited_in_walk.contains(&e.to))
            .collect();

        let chosen: Edge;
        let is_revisit: bool;
        if !unvisited.is_empty() {
            // Pick among unvisited, weighted to prefer non-global-avoid neighbors.
            let weights: Vec<f64> = unvisited
                .iter()
                .map(|e| {
                    let mut w = 1.0;
                    if let Some(av) = global_avoid {
                        if av.contains(&e.to) {
                            w /= global_avoid_weight.max(1.0);
                        }
                    }
                    w
                })
                .collect();
            let idx = rng.weighted_pick(&weights);
            chosen = unvisited[idx].clone();
            is_revisit = false;
        } else {
            // No unvisited neighbors. Are we allowed to revisit?
            if revisits_used >= allow_revisits {
                break;
            }
            // Pick any edge (still bias against global_avoid).
            let weights: Vec<f64> = edges
                .iter()
                .map(|e| {
                    let mut w = 1.0;
                    if let Some(av) = global_avoid {
                        if av.contains(&e.to) {
                            w /= global_avoid_weight.max(1.0);
                        }
                    }
                    w
                })
                .collect();
            let idx = rng.weighted_pick(&weights);
            chosen = edges[idx].clone();
            is_revisit = true;
        }

        let next_node = chosen.to;
        let next_entry = chosen.to_side;
        let next_exit: Side = 1 - next_entry;
        let next_orient = if next_entry == 0 { '+' } else { '-' };

        walk.steps.push(Step {
            node: next_node,
            orient: next_orient,
            overlap_in: chosen.overlap,
            cigar_in: chosen.cigar.clone(),
        });
        if is_revisit {
            revisits_used += 1;
        }
        visited_in_walk.insert(next_node);

        cur = next_node;
        cur_exit = next_exit;
    }

    walk
}

/// Generate K candidate walks for one component using several heuristic
/// flavors and seeds. The pool is then deduplicated lightly by node-set
/// identity.
fn generate_candidates(
    g: &Graph,
    members: &[usize],
    mode: EndpointMode,
    k: usize,
    base_seed: u64,
) -> Vec<Walk> {
    // Determine starts under the requested mode, with tip->any fallback.
    let mut starts = walk_starts(g, members, mode);
    if starts.is_empty() {
        starts = walk_starts(g, members, EndpointMode::Any);
    }
    if starts.is_empty() {
        return Vec::new();
    }

    let mut pool: Vec<Walk> = Vec::with_capacity(k);

    // Mix of revisit budgets: try strict simple-path walks (0 revisits) AND
    // walks with small budgets that can dip back to grab extra coverage.
    let revisit_budgets: [usize; 4] = [0, 1, 3, 8];

    let mut seed_counter: u64 = 0;
    for trial in 0..k {
        let budget = revisit_budgets[trial % revisit_budgets.len()];
        let start_idx = trial % starts.len();
        let (sn, se) = starts[start_idx];

        let mut rng = Rng::new(
            base_seed
                .wrapping_add(seed_counter)
                .wrapping_mul(0x100000001B3),
        );
        seed_counter = seed_counter.wrapping_add(1);

        let walk = generate_walk(g, sn, se, &mut rng, budget, None, 1.0);
        pool.push(walk);
    }

    // Dedup by node-multiset signature.
    let mut seen_sigs: HashSet<Vec<usize>> = HashSet::new();
    let mut deduped: Vec<Walk> = Vec::with_capacity(pool.len());
    for w in pool {
        let mut sig: Vec<usize> = w.steps.iter().map(|s| s.node).collect();
        sig.sort_unstable();
        if seen_sigs.insert(sig) {
            deduped.push(w);
        }
    }
    deduped
}

/// Given a pool of candidates, find the best pair by lexicographic objective.
/// Score returned alongside for reporting.
#[derive(Clone, Copy, Debug)]
struct PairScore {
    union_bp: u64,
    symdiff_bp: u64,
    revisits_total: i64, // negated for minimization in lex comparison
}

fn pair_score(g: &Graph, w1: &Walk, w2: &Walk) -> PairScore {
    let s1 = w1.node_set();
    let s2 = w2.node_set();
    let mut union_bp: u64 = 0;
    let mut sym_bp: u64 = 0;
    let union: HashSet<usize> = s1.union(&s2).copied().collect();
    for &n in &union {
        let len = g.segs[n].length;
        union_bp += len;
        let in1 = s1.contains(&n);
        let in2 = s2.contains(&n);
        if in1 ^ in2 {
            sym_bp += len;
        }
    }
    let r = (w1.revisits() + w2.revisits()) as i64;
    PairScore {
        union_bp,
        symdiff_bp: sym_bp,
        revisits_total: r,
    }
}

fn pair_better(a: PairScore, b: PairScore) -> bool {
    // Lex: bigger union, then bigger symdiff, then fewer revisits.
    if a.union_bp != b.union_bp {
        return a.union_bp > b.union_bp;
    }
    if a.symdiff_bp != b.symdiff_bp {
        return a.symdiff_bp > b.symdiff_bp;
    }
    a.revisits_total < b.revisits_total
}

/// Top-level: generate candidates, optionally augment with "diverging" walks
/// computed against the best-so-far candidate, then pick the best pair.
fn best_pair_for_component(
    g: &Graph,
    members: &[usize],
    mode: EndpointMode,
    k: usize,
    base_seed: u64,
    verbose: bool,
) -> Option<(Walk, Walk, PairScore)> {
    let mut pool = generate_candidates(g, members, mode, k, base_seed);
    if pool.is_empty() {
        return None;
    }

    // Find the best single walk by length first; use its node set to guide
    // generation of "complementary" walks.
    pool.sort_by(|a, b| b.length_bp(g).cmp(&a.length_bp(g)));
    let leader = pool[0].clone();
    let leader_nodes: HashSet<usize> = leader.node_set();

    // Generate complementary walks that avoid leader's nodes (haplotype-divergence guide).
    let mut starts = walk_starts(g, members, mode);
    if starts.is_empty() {
        starts = walk_starts(g, members, EndpointMode::Any);
    }
    let mut comp_seed: u64 = base_seed ^ 0xA5A5_5A5A_DEAD_BEEF;
    let comp_count = k; // same budget for complements
    let revisit_budgets: [usize; 4] = [0, 1, 3, 8];
    let global_avoid_weights: [f64; 4] = [4.0, 16.0, 64.0, 256.0];
    for trial in 0..comp_count {
        let budget = revisit_budgets[trial % revisit_budgets.len()];
        let gw = global_avoid_weights[trial % global_avoid_weights.len()];
        let (sn, se) = starts[trial % starts.len()];
        let mut rng = Rng::new(comp_seed);
        comp_seed = comp_seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let w = generate_walk(g, sn, se, &mut rng, budget, Some(&leader_nodes), gw);
        pool.push(w);
    }

    // Evaluate all pairs (ordered as the canonical (smaller_idx, larger_idx)).
    if verbose {
        eprintln!("  pool size: {}", pool.len());
        // Print top candidates by length
        let mut idx: Vec<usize> = (0..pool.len()).collect();
        idx.sort_by_key(|&i| std::cmp::Reverse(pool[i].length_bp(g)));
        for &i in idx.iter().take(20) {
            let ns: HashSet<usize> = pool[i].node_set();
            eprintln!(
                "    cand[{}] len_bp={} nodes={} revisits={} walk={}",
                i,
                pool[i].length_bp(g),
                ns.len(),
                pool[i].revisits(),
                pool[i].to_walk_string(g)
            );
        }
    }
    let mut best: Option<(usize, usize, PairScore)> = None;
    for i in 0..pool.len() {
        for j in (i + 1)..pool.len() {
            let sc = pair_score(g, &pool[i], &pool[j]);
            match best {
                None => best = Some((i, j, sc)),
                Some((_, _, bs)) => {
                    if pair_better(sc, bs) {
                        best = Some((i, j, sc));
                    }
                }
            }
        }
        // Also allow degenerate "pair with itself" only if pool has size 1
        if pool.len() == 1 {
            let sc = pair_score(g, &pool[i], &pool[i]);
            match best {
                None => best = Some((i, i, sc)),
                Some(_) => {}
            }
        }
    }

    best.map(|(i, j, sc)| (pool[i].clone(), pool[j].clone(), sc))
}

// ------------ Coverage repair ------------

/// Validate that a walk is structurally well-formed:
/// 1. No adjacent identical steps (same node, same orientation).
/// 2. Each consecutive (step_i, step_{i+1}) pair is connected by a real
///    edge: exiting step_i on its exit-side leads to entering step_{i+1}
///    on its entry-side via some L-line in the graph.
fn is_valid_walk(g: &Graph, walk: &Walk) -> bool {
    if walk.steps.is_empty() {
        return true;
    }
    for i in 1..walk.steps.len() {
        let prev = &walk.steps[i - 1];
        let cur = &walk.steps[i];
        // Reject adjacent same-node same-orientation (the bug we hit).
        // (Different orientations of the same node back-to-back are still
        // suspicious in a unitig walk, but they could in principle arise
        // from a hairpin; we don't reject those, only true duplicates.)
        if prev.node == cur.node && prev.orient == cur.orient {
            return false;
        }
        let (_, prev_exit) = entry_exit_for_orient(prev.orient);
        let (cur_entry, _) = entry_exit_for_orient(cur.orient);
        // There must be an edge from (prev.node, prev_exit) to
        // (cur.node, cur_entry).
        let found = g
            .neighbors(prev.node, prev_exit)
            .iter()
            .any(|e| e.to == cur.node && e.to_side == cur_entry);
        if !found {
            return false;
        }
    }
    true
}

// After best-pair selection, attempt to extend coverage by splicing detours
// for uncovered nodes into one of the two walks. A detour through node X is
// a path:  walk_step_i  -> ... -> X -> ... -> walk_step_j  (j >= i+1)
// that we insert in place of the (i, i+1) edge. The replacement must:
//   * be a valid bidirected path in the graph
//   * leave walk_step_i exiting from the same side it would have under the
//     original walk (otherwise the rest of the walk is broken)
//   * arrive at walk_step_j entering from the same side it originally did
//   * (we DO allow re-using nodes already in the walk, but we count those
//     as revisits.)
//
// Strategy:
// - Build a BFS from each (node, exit_side) state.
// - For each uncovered node U and each side S of U we can enter, search for
//   a state in the "destinations" reachable from U that matches some
//   subsequent step of the walk. If found and the entry-side into U is
//   reachable from the corresponding walk's predecessor exit-state, splice.
// - Try both walks; prefer the splice with fewer added revisits.
//
// This is necessarily heuristic for large graphs (full bidirected
// reachability with side states is up to 2N states, BFS is O(E) each).
// We cap effort per uncovered node.

use std::collections::VecDeque;

/// BFS over the bidirected state graph (node, exit_side). Returns parent
/// pointers usable to reconstruct paths. Source state is (start_node,
/// start_exit_side): "I have just entered start_node and am about to exit
/// via start_exit_side". The BFS expands: take an outgoing edge to (to,
/// to_entry_side); the next state is (to, 1 - to_entry_side) -- i.e., we
/// then exit `to` via its opposite side.
///
/// Returns parents indexed by halfedge(node, exit_side), each entry
/// (parent_halfedge, overlap_to_get_here, cigar_to_get_here) or None.
fn bfs_from(
    g: &Graph,
    start_node: usize,
    start_exit_side: Side,
) -> Vec<Option<(usize, u64, String)>> {
    let n = g.n();
    let mut parent: Vec<Option<(usize, u64, String)>> = vec![None; n * 2];
    let mut q: VecDeque<usize> = VecDeque::new();
    let s = halfedge(start_node, start_exit_side);
    // Sentinel for source: parent index = usize::MAX.
    parent[s] = Some((usize::MAX, 0, String::new()));
    q.push_back(s);
    while let Some(h) = q.pop_front() {
        let node = h / 2;
        let exit_side = (h % 2) as Side;
        for e in g.neighbors(node, exit_side) {
            let next_exit = 1 - e.to_side;
            let nh = halfedge(e.to, next_exit);
            if parent[nh].is_none() {
                parent[nh] = Some((h, e.overlap, e.cigar.clone()));
                q.push_back(nh);
            }
        }
    }
    parent
}

/// Reconstruct path from BFS parents. Each entry is
/// (node, exit_side, overlap_used_to_arrive, cigar_used_to_arrive).
/// The first entry is the source state with overlap=0, cigar="".
fn reconstruct_path(
    parents: &[Option<(usize, u64, String)>],
    target_he: usize,
) -> Option<Vec<(usize, Side, u64, String)>> {
    parents.get(target_he)?.as_ref()?;
    let mut path: Vec<(usize, Side, u64, String)> = Vec::new();
    let mut cur = target_he;
    loop {
        let entry = parents[cur].as_ref().unwrap();
        let par = entry.0;
        let ov = entry.1;
        let cig = entry.2.clone();
        let node = cur / 2;
        let exit_side = (cur % 2) as Side;
        path.push((node, exit_side, ov, cig));
        if par == usize::MAX {
            break;
        }
        cur = par;
    }
    path.reverse();
    Some(path)
}

/// Convert a path-of-exit-states into Steps for splicing into a Walk.
/// path[0] is the splice predecessor and is NOT emitted; emitted steps are
/// path[1..].
fn path_to_steps(g: &Graph, path: &[(usize, Side, u64, String)]) -> Vec<Step> {
    let mut out = Vec::new();
    for i in 1..path.len() {
        let (this_node, this_exit, ov, ref cig) = path[i];
        let entry_side = 1 - this_exit;
        let orient = if entry_side == 0 { '+' } else { '-' };
        let _ = g;
        out.push(Step {
            node: this_node,
            orient,
            overlap_in: ov,
            cigar_in: cig.clone(),
        });
    }
    out
}

/// Build the "exit state" after step `i` in the walk: the (node, exit_side)
/// from which the walk continues to step i+1.
fn exit_state_after(walk: &Walk, i: usize) -> (usize, Side) {
    let step = &walk.steps[i];
    let (_, exit) = entry_exit_for_orient(step.orient);
    (step.node, exit)
}

/// Build the "entry state" expected for step `i`: the (node, entry_side)
/// the predecessor must hand us so that this step's orientation is preserved.
fn entry_state_of(walk: &Walk, i: usize) -> (usize, Side) {
    let step = &walk.steps[i];
    let (entry, _) = entry_exit_for_orient(step.orient);
    (step.node, entry)
}

/// Attempt to splice node `target` into `walk` such that walk visits target
/// somewhere. Returns Some(new_walk) on success, None if no splice found.
///
/// We try every (i, j) split where 0 <= i < j <= len-1, finding a path from
/// exit-state-after(i) through `target` to entry-state-of(j). If multiple
/// splices succeed, return the one with the fewest *new revisits* (nodes
/// added to the walk that were already present in the walk).
fn try_splice(g: &Graph, walk: &Walk, target: usize, max_extra_revisits: usize) -> Option<Walk> {
    let n_steps = walk.steps.len();
    if n_steps < 1 {
        return None;
    }

    // Pre-build BFS from each step's exit state (memoize across i's).
    let mut bfs_cache: HashMap<usize, Vec<Option<(usize, u64, String)>>> = HashMap::new();

    let original_nodes: HashSet<usize> = walk.steps.iter().map(|s| s.node).collect();

    let mut best: Option<(Walk, usize)> = None;

    // Helper that gets or computes a BFS result, returning a clone of the vec.
    // (Cloning is acceptable here: per call site we do small amounts of work
    // on the result and the per-component BFS is cheap relative to candidate
    // generation.)
    fn get_bfs(
        g: &Graph,
        cache: &mut HashMap<usize, Vec<Option<(usize, u64, String)>>>,
        node: usize,
        exit_side: Side,
    ) -> Vec<Option<(usize, u64, String)>> {
        let key = halfedge(node, exit_side);
        if let Some(v) = cache.get(&key) {
            return v.clone();
        }
        let v = bfs_from(g, node, exit_side);
        cache.insert(key, v.clone());
        v
    }

    for i in 0..n_steps {
        let (a_node, a_exit) = exit_state_after(walk, i);
        let parents_fwd = get_bfs(g, &mut bfs_cache, a_node, a_exit);

        for target_entry in [0u8, 1u8] {
            let target_exit_state = halfedge(target, 1 - target_entry);
            if parents_fwd[target_exit_state].is_none() {
                continue;
            }
            let path_to_target = match reconstruct_path(&parents_fwd, target_exit_state) {
                Some(p) => p,
                None => continue,
            };
            let parents_from_t = get_bfs(g, &mut bfs_cache, target, 1 - target_entry);

            for j in (i + 1)..n_steps {
                let (b_node, b_entry) = entry_state_of(walk, j);
                let target_he = halfedge(b_node, 1 - b_entry);
                if parents_from_t[target_he].is_none() {
                    continue;
                }
                let path_from_target = match reconstruct_path(&parents_from_t, target_he) {
                    Some(p) => p,
                    None => continue,
                };

                let mut detour_combined: Vec<(usize, Side, u64, String)> =
                    path_to_target.clone();
                detour_combined.extend(path_from_target.iter().skip(1).cloned());

                let detour_steps = path_to_steps(g, &detour_combined);

                // The last entry of detour_combined is the state (b_node,
                // 1 - b_entry) -- i.e., "we just entered b_node on b_entry,
                // about to exit". That state IS walk.steps[j]: same node,
                // same entry side, so same orientation. We must therefore
                // skip walk.steps[j] when rebuilding, otherwise it appears
                // twice (an adjacent duplicate, which breaks Bandage and
                // produces an invalid walk).
                //
                // We also need to preserve walk.steps[j]'s overlap_in for the
                // detour's last step, since that overlap is what bridges to
                // the splice point. Actually no -- the detour's last step
                // already carries the correct overlap from path_from_target's
                // final edge, which is the edge we actually traverse in the
                // spliced walk. walk.steps[j].overlap_in was the overlap in
                // the *original* walk (from j-1 to j), which no longer
                // applies after splicing.

                let mut seen: HashSet<usize> = original_nodes.clone();
                let mut extras = 0usize;
                for st in &detour_steps {
                    if !seen.insert(st.node) {
                        extras += 1;
                    }
                }

                if extras > max_extra_revisits {
                    continue;
                }

                let mut new_steps: Vec<Step> = walk.steps[..=i].to_vec();
                new_steps.extend(detour_steps);
                // Skip walk.steps[j] -- it's already the last detour step.
                if j + 1 < walk.steps.len() {
                    new_steps.extend(walk.steps[(j + 1)..].iter().cloned());
                }

                let new_walk = Walk { steps: new_steps };

                let included = new_walk.steps.iter().any(|s| s.node == target);
                if !included {
                    continue;
                }

                // Safety check: reject any walk that produces adjacent
                // duplicate nodes or whose adjacent steps aren't connected
                // by a real edge. This guards against any remaining splice
                // edge-cases we haven't enumerated.
                if !is_valid_walk(g, &new_walk) {
                    continue;
                }

                match &best {
                    None => best = Some((new_walk, extras)),
                    Some((_, prev_extras)) => {
                        if extras < *prev_extras {
                            best = Some((new_walk, extras));
                        }
                    }
                }
            }
        }
    }

    best.map(|(w, _)| w)
}

/// Try to extend a walk past its current end to pick up an uncovered node.
/// This is a simpler variant of try_splice that only appends to the end of
/// the walk. Useful when the walk currently terminates at an internal node
/// because the algorithm decided to stop.
fn try_extend_end(g: &Graph, walk: &Walk, target: usize, max_extra_revisits: usize) -> Option<Walk> {
    let n_steps = walk.steps.len();
    if n_steps < 1 {
        return None;
    }
    let last = &walk.steps[n_steps - 1];
    let (_, exit) = entry_exit_for_orient(last.orient);
    let parents = bfs_from(g, last.node, exit);
    let original_nodes: HashSet<usize> = walk.steps.iter().map(|s| s.node).collect();

    let mut best: Option<(Walk, usize)> = None;
    for target_entry in [0u8, 1u8] {
        let target_exit_state = halfedge(target, 1 - target_entry);
        if parents[target_exit_state].is_none() {
            continue;
        }
        let path = match reconstruct_path(&parents, target_exit_state) {
            Some(p) => p,
            None => continue,
        };
        let append_steps = path_to_steps(g, &path);
        let mut seen = original_nodes.clone();
        let mut extras = 0;
        for st in &append_steps {
            if !seen.insert(st.node) {
                extras += 1;
            }
        }
        if extras > max_extra_revisits {
            continue;
        }
        let mut new_steps = walk.steps.clone();
        new_steps.extend(append_steps);
        let new_walk = Walk { steps: new_steps };
        if !is_valid_walk(g, &new_walk) {
            continue;
        }
        match &best {
            None => best = Some((new_walk, extras)),
            Some((_, pe)) => {
                if extras < *pe {
                    best = Some((new_walk, extras));
                }
            }
        }
    }
    best.map(|(w, _)| w)
}

/// Repair coverage for a pair of walks. Tries to add each uncovered node
/// into one of the two walks via try_splice or try_extend_end, choosing the
/// option that yields the highest *joint* coverage afterward. Splices are
/// only accepted if they strictly increase joint coverage, preventing
/// flip-flops where splicing one node evicts another.
fn repair_coverage(
    g: &Graph,
    members: &[usize],
    w1: Walk,
    w2: Walk,
    max_extra_revisits_per_node: usize,
    verbose: bool,
) -> (Walk, Walk, usize, usize) {
    let mut w1 = w1;
    let mut w2 = w2;
    let mut splices = 0usize;
    // Hard iteration cap to guarantee termination even if something pathological happens.
    let mut iters = 0usize;
    let max_iters = members.len() * 4 + 16;

    loop {
        iters += 1;
        if iters > max_iters {
            if verbose {
                eprintln!("    repair hit iteration cap ({})", max_iters);
            }
            break;
        }

        let covered_before: HashSet<usize> =
            w1.node_set().union(&w2.node_set()).copied().collect();
        let uncovered: Vec<usize> = members
            .iter()
            .copied()
            .filter(|u| !covered_before.contains(u))
            .collect();
        if uncovered.is_empty() {
            break;
        }

        let cov_before_count = covered_before.len();

        // For each uncovered node, find the best splice option (if any) and
        // remember the resulting (which_walk, new_walk, new_joint_coverage).
        // Apply the single best splice this iteration.
        let mut best: Option<(usize, usize, Walk, usize)> = None;
        // (target_node, which_walk_index (0 or 1), new_walk, new_coverage_count)

        for u in &uncovered {
            for which in 0..2usize {
                let w_ref = if which == 0 { &w1 } else { &w2 };
                let new_w_opt = try_splice(g, w_ref, *u, max_extra_revisits_per_node)
                    .or_else(|| try_extend_end(g, w_ref, *u, max_extra_revisits_per_node));
                if let Some(new_w) = new_w_opt {
                    // Compute joint coverage if we replace which-th walk with new_w
                    let other = if which == 0 { &w2 } else { &w1 };
                    let joint: HashSet<usize> =
                        new_w.node_set().union(&other.node_set()).copied().collect();
                    let joint_n = joint.len();
                    // Only consider if joint coverage strictly improves.
                    if joint_n > cov_before_count {
                        let should_take = match &best {
                            None => true,
                            Some((_, _, _, best_n)) => joint_n > *best_n,
                        };
                        if should_take {
                            best = Some((*u, which, new_w, joint_n));
                        }
                    }
                }
            }
        }

        match best {
            Some((target, which, new_w, joint_n)) => {
                if verbose {
                    eprintln!(
                        "    splice {} into w{} (coverage {} -> {})",
                        g.segs[target].name,
                        which + 1,
                        cov_before_count,
                        joint_n
                    );
                }
                if which == 0 {
                    w1 = new_w;
                } else {
                    w2 = new_w;
                }
                splices += 1;
            }
            None => {
                if verbose {
                    eprintln!(
                        "    coverage repair stalled with {} node(s) still uncovered",
                        uncovered.len()
                    );
                }
                break;
            }
        }
    }

    let final_uncovered = members
        .iter()
        .filter(|u| {
            !w1.node_set().contains(u) && !w2.node_set().contains(u)
        })
        .count();
    (w1, w2, splices, final_uncovered)
}



fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'N' | b'n' => b,
            _ => b'N',
        })
        .collect()
}

/// Stitch a walk into a sequence, trimming overlaps off the prefix of each
/// non-first segment. Returns None if any segment lacks a stored sequence.
fn walk_to_sequence(g: &Graph, walk: &Walk) -> Option<Vec<u8>> {
    let mut out: Vec<u8> = Vec::new();
    for (i, st) in walk.steps.iter().enumerate() {
        let seg = &g.segs[st.node];
        let seq = seg.seq.as_ref()?;
        let oriented: Vec<u8> = if st.orient == '+' {
            seq.clone()
        } else {
            revcomp(seq)
        };
        if i == 0 {
            out.extend_from_slice(&oriented);
        } else {
            let trim = st.overlap_in as usize;
            if trim >= oriented.len() {
                // pathological: overlap >= segment length; skip
                continue;
            }
            out.extend_from_slice(&oriented[trim..]);
        }
    }
    Some(out)
}

// ------------ Main ------------

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "usage: {} <input.gfa> [--endpoints tip|any] [--candidates K] [--seed S] [--sample NAME] [--skip-sequences] [--repair-max-revisits N] [--no-repair] [--verbose]",
            args.get(0).map(|s| s.as_str()).unwrap_or("gfa_haps")
        );
        process::exit(2);
    }

    let input_path = PathBuf::from(&args[1]);
    let mut mode = EndpointMode::Tip;
    let mut k: usize = 64;
    let mut seed: u64 = 0xC0FFEE;
    let mut sample = String::from("sample");
    let mut skip_sequences = false;
    let mut verbose = false;
    let mut repair_max_revisits: usize = 4;
    let mut no_repair = false;

    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--endpoints" => {
                i += 1;
                if i >= args.len() {
                    die("--endpoints needs a value");
                }
                mode = match args[i].as_str() {
                    "tip" => EndpointMode::Tip,
                    "any" => EndpointMode::Any,
                    other => die(&format!("--endpoints: unknown value {}", other)),
                };
            }
            "--candidates" => {
                i += 1;
                if i >= args.len() {
                    die("--candidates needs a value");
                }
                k = args[i].parse().unwrap_or_else(|_| die("--candidates must be a positive integer"));
                if k < 2 {
                    die("--candidates must be >= 2");
                }
            }
            "--seed" => {
                i += 1;
                if i >= args.len() {
                    die("--seed needs a value");
                }
                seed = args[i].parse().unwrap_or_else(|_| die("--seed must be an integer"));
            }
            "--sample" => {
                i += 1;
                if i >= args.len() {
                    die("--sample needs a value");
                }
                sample = args[i].clone();
            }
            "--skip-sequences" => {
                skip_sequences = true;
            }
            "--verbose" | "-v" => {
                verbose = true;
            }
            "--repair-max-revisits" => {
                i += 1;
                if i >= args.len() {
                    die("--repair-max-revisits needs a value");
                }
                repair_max_revisits = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--repair-max-revisits must be a non-negative integer"));
            }
            "--no-repair" => {
                no_repair = true;
            }
            other => die(&format!("unknown argument: {}", other)),
        }
        i += 1;
    }

    let opts = ParseOptions {
        keep_sequences: true,
    };

    eprintln!("parsing {}", input_path.display());
    let parsed = parse_gfa(&input_path, &opts).unwrap_or_else(|e| die(&format!("parse error: {}", e)));
    let g = &parsed.graph;
    eprintln!("  {} segments", g.n());

    let (comp_of, n_comp, members) = connected_components(g);
    let _ = comp_of;
    eprintln!("  {} connected component(s)", n_comp);

    // Sort components by total seg length descending so output order is stable
    // and comp_01 is the biggest.
    let mut order: Vec<usize> = (0..n_comp).collect();
    order.sort_by_key(|&ci| {
        let total: u64 = members[ci].iter().map(|&u| g.segs[u].length).sum();
        std::cmp::Reverse(total)
    });

    // Open outputs
    let stem = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("out")
        .to_string();
    let out_dir = input_path
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));
    let gfa_out_path = out_dir.join(format!("{}.with_paths.gfa", stem));
    let fa_out_path = out_dir.join(format!("{}.haps.fa", stem));

    let mut gfa_out =
        BufWriter::new(File::create(&gfa_out_path).unwrap_or_else(|e| {
            die(&format!("cannot create {}: {}", gfa_out_path.display(), e))
        }));
    let mut fa_out = BufWriter::new(File::create(&fa_out_path).unwrap_or_else(
        |e| die(&format!("cannot create {}: {}", fa_out_path.display(), e)),
    ));

    // First, re-stream the input to copy all original lines (minus A-lines)
    // to the output. We don't try to merge with the W-lines we'll add; we
    // simply append our new W-lines at the end.
    let f = File::open(&input_path).unwrap_or_else(|e| die(&format!("cannot reopen input: {}", e)));
    use std::io::{BufRead, BufReader};
    let reader = BufReader::new(f);
    for line in reader.lines() {
        let line = line.unwrap_or_else(|e| die(&format!("read error: {}", e)));
        if line.is_empty() {
            continue;
        }
        // Skip A-lines on output, same convention as gfa_split.
        let rec_end = line.find('\t').unwrap_or(line.len());
        let rec = &line[..rec_end];
        if rec == "A" {
            continue;
        }
        if skip_sequences && rec == "S" {
            // Replace the sequence field with "*". Split into at most 4 parts
            // so that tags (field 3+) are kept as-is. If the sequence was
            // non-"*", add an LN:i: tag so consumers can still recover the
            // segment length.
            let fields: Vec<&str> = line.splitn(4, '\t').collect();
            if fields.len() >= 3 && fields[2] != "*" {
                let seq_len = fields[2].len();
                let tags = fields.get(3).copied().unwrap_or("");
                let has_ln = tags.split('\t').any(|t| t.starts_with("LN:i:"));
                if tags.is_empty() {
                    writeln!(&mut gfa_out, "S\t{}\t*\tLN:i:{}", fields[1], seq_len).unwrap();
                } else if has_ln {
                    writeln!(&mut gfa_out, "S\t{}\t*\t{}", fields[1], tags).unwrap();
                } else {
                    writeln!(&mut gfa_out, "S\t{}\t*\t{}\tLN:i:{}", fields[1], tags, seq_len).unwrap();
                }
                continue;
            }
        }
        writeln!(&mut gfa_out, "{}", line).unwrap();
    }

    // Per-component summary
    let pad = std::cmp::max(2, n_comp.to_string().len());
    println!("# component\tsize_bp\thap1_len\thap2_len\tunion_bp\tsymdiff_bp\trevisits\tnodes_total\tnodes_covered");
    for (out_idx0, &ci) in order.iter().enumerate() {
        let label_idx = out_idx0 + 1;
        let m = &members[ci];
        let component_label = format!("component_{:0width$}", label_idx, width = pad);
        let total_bp: u64 = m.iter().map(|&u| g.segs[u].length).sum();

        let comp_seed = seed.wrapping_add((ci as u64).wrapping_mul(0x9E3779B97F4A7C15));
        match best_pair_for_component(g, m, mode, k, comp_seed, verbose) {
            None => {
                eprintln!("warning: no walks generated for {}", component_label);
            }
            Some((w1, w2, _sc)) => {
                // Coverage repair: try to splice uncovered nodes into walks.
                let (w1, w2, splices, still_uncovered) = if no_repair {
                    (w1, w2, 0usize, 0usize)
                } else {
                    repair_coverage(g, m, w1, w2, repair_max_revisits, verbose)
                };
                if verbose && splices > 0 {
                    eprintln!(
                        "  repair: {} splice(s), {} node(s) still uncovered",
                        splices, still_uncovered
                    );
                }
                // Recompute score after possible repair
                let sc = pair_score(g, &w1, &w2);
                let s1 = w1.node_set();
                let s2 = w2.node_set();
                let covered: HashSet<usize> = s1.union(&s2).copied().collect();
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    component_label,
                    total_bp,
                    w1.length_bp(g),
                    w2.length_bp(g),
                    sc.union_bp,
                    sc.symdiff_bp,
                    sc.revisits_total,
                    m.len(),
                    covered.len()
                );

                // Emit P-lines (Path):
                //   P <path_name> <seg1+,seg2+,...> <cigar1,cigar2,...>
                //
                // We use P-lines (not W-lines) because the GFA spec
                // restricts W-lines to graphs without overlaps between
                // segments. Hifiasm unitig graphs have nonzero overlaps on
                // every link, so P-lines are the correct format here.
                //
                // Path name encodes sample, haplotype, and component label
                // so two haplotype paths are distinguishable.
                //
                // Final validation -- this should always pass, but if a bug
                // ever produces a malformed walk, we'd rather fail loudly
                // than write a broken GFA file.
                if !is_valid_walk(g, &w1) {
                    eprintln!(
                        "error: {} hap1 is structurally invalid; this is a bug. Skipping.",
                        component_label
                    );
                    continue;
                }
                if !is_valid_walk(g, &w2) {
                    eprintln!(
                        "error: {} hap2 is structurally invalid; this is a bug. Skipping.",
                        component_label
                    );
                    continue;
                }
                let segs1 = w1.to_path_segment_names(g);
                let ovs1 = w1.to_path_overlaps();
                let segs2 = w2.to_path_segment_names(g);
                let ovs2 = w2.to_path_overlaps();
                writeln!(
                    &mut gfa_out,
                    "P\t{}_h1_{}\t{}\t{}",
                    sample, component_label, segs1, ovs1
                )
                .unwrap();
                writeln!(
                    &mut gfa_out,
                    "P\t{}_h2_{}\t{}\t{}",
                    sample, component_label, segs2, ovs2
                )
                .unwrap();

                if let Some(seq) = walk_to_sequence(g, &w1) {
                    writeln!(fa_out, ">{}_hap1 len={}", component_label, seq.len()).unwrap();
                    write_wrapped(&mut fa_out, &seq, 80);
                }
                if let Some(seq) = walk_to_sequence(g, &w2) {
                    writeln!(fa_out, ">{}_hap2 len={}", component_label, seq.len()).unwrap();
                    write_wrapped(&mut fa_out, &seq, 80);
                }
            }
        }
    }

    gfa_out.flush().unwrap();
    fa_out.flush().unwrap();
    eprintln!("wrote {} and {}", gfa_out_path.display(), fa_out_path.display());
}

fn write_wrapped<W: Write>(w: &mut W, seq: &[u8], width: usize) {
    let mut i = 0;
    while i < seq.len() {
        let end = (i + width).min(seq.len());
        w.write_all(&seq[i..end]).unwrap();
        w.write_all(b"\n").unwrap();
        i = end;
    }
}
