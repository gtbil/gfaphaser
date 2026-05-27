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

mod diverse_pair;
mod divide_and_conquer;

use diverse_pair::{pair_score, DiversePair, EndpointMode, PairScore};
use divide_and_conquer::{collect_subgraph_nodes, find_anchors, order_anchors, DivideAndConquer};
use rayon;
use rayon::prelude::*;
use gfaphaser::{
    connected_components, entry_exit_for_orient, halfedge, parse_gfa, Graph, ParseOptions,
    Side, Step, Walk,
};
use std::collections::{HashMap, HashSet, VecDeque};
use std::rc::Rc;

// BFS parent table keyed by halfedge(node, exit_side).
// Values are Rc so "cloning" a cached entry is O(1) (just a refcount bump),
// not a deep copy of the O(N) parent vector.
type BfsEntry = Rc<Vec<Option<(usize, u64, String)>>>;
type BfsCache = HashMap<usize, BfsEntry>;

fn get_bfs(g: &Graph, cache: &mut BfsCache, node: usize, exit_side: Side) -> BfsEntry {
    let key = halfedge(node, exit_side);
    if let Some(v) = cache.get(&key) {
        return Rc::clone(v);
    }
    let v = Rc::new(bfs_from(g, node, exit_side));
    cache.insert(key, Rc::clone(&v));
    v
}
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process;

fn die(msg: &str) -> ! {
    eprintln!("error: {}", msg);
    process::exit(1);
}

// ------------ Walk-pair algorithm trait ------------

trait WalkPairAlgorithm {
    fn build_pair(
        &self,
        g: &Graph,
        members: &[usize],
        comp_idx: usize, // original (pre-sort) component index; used for per-component seed derivation
        verbose: bool,
    ) -> Option<(Walk, Walk)>;
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
/// Searches in order of increasing gap (j - i) so that short detours are
/// found first. A short splice replaces few existing steps, preserving the
/// overall walk structure and reducing the risk of evicting nodes that are
/// only covered by this walk. We return as soon as the first valid splice
/// is found — no need to scan all pairs.
///
/// Complexity: O(G × n_steps) BFS lookups (O(1) each from cache) where G
/// is the gap of the first valid splice found. Path reconstruction and walk
/// building are O(path_len) and happen only for the winning pair.
fn try_splice(
    g: &Graph,
    walk: &Walk,
    target: usize,
    max_extra_revisits: usize,
    bfs_cache: &mut BfsCache,
) -> Option<Walk> {
    let n_steps = walk.steps.len();
    if n_steps < 2 {
        return None;
    }

    let original_nodes: HashSet<usize> = walk.steps.iter().map(|s| s.node).collect();

    // Precompute: for each target_entry side (0 or 1), the set of walk step
    // indices j whose entry state is reachable from target's exit. Stored as
    // HashSets for O(1) per-j filter inside the gap loop.
    let reachable_j: [HashSet<usize>; 2] = [0u8, 1u8].map(|te| {
        let parents_from_t = get_bfs(g, bfs_cache, target, 1 - te);
        (0..n_steps)
            .filter(|&j| {
                let (b_node, b_entry) = entry_state_of(walk, j);
                parents_from_t[halfedge(b_node, 1 - b_entry)].is_some()
            })
            .collect::<HashSet<usize>>()
    });

    // If target can reach no step j at all, bail immediately.
    if reachable_j[0].is_empty() && reachable_j[1].is_empty() {
        return None;
    }

    // Scan increasing gap values. For gap g, try all (i, j = i+g) pairs.
    // Short splices are less disruptive so we prefer them by trying first.
    for gap in 1..n_steps {
        for i in 0..(n_steps - gap) {
            let j = i + gap;
            for target_entry in [0u8, 1u8] {
                // O(1) filter: j reachable from target?
                if !reachable_j[target_entry as usize].contains(&j) {
                    continue;
                }
                // O(1) filter: target reachable from walk[i]'s exit?
                let (a_node, a_exit) = exit_state_after(walk, i);
                let parents_fwd = get_bfs(g, bfs_cache, a_node, a_exit);
                let target_exit_state = halfedge(target, 1 - target_entry);
                if parents_fwd[target_exit_state].is_none() {
                    continue;
                }

                // Both reachable — reconstruct paths and check revisit budget.
                let path_to_target = match reconstruct_path(&parents_fwd, target_exit_state) {
                    Some(p) => p,
                    None => continue,
                };
                let parents_from_t = get_bfs(g, bfs_cache, target, 1 - target_entry);
                let (b_node, b_entry) = entry_state_of(walk, j);
                let target_he = halfedge(b_node, 1 - b_entry);
                let path_from_target = match reconstruct_path(&parents_from_t, target_he) {
                    Some(p) => p,
                    None => continue,
                };

                let mut seen: HashSet<usize> = HashSet::new();
                let mut extras = 0usize;
                for &(n, _, _, _) in path_to_target[1..]
                    .iter()
                    .chain(path_from_target[1..].iter())
                {
                    if original_nodes.contains(&n) || !seen.insert(n) {
                        extras += 1;
                    }
                }
                if extras > max_extra_revisits {
                    continue;
                }

                // Build the detour and splice it in.
                let mut detour: Vec<(usize, Side, u64, String)> = path_to_target;
                detour.extend(path_from_target.into_iter().skip(1));
                let detour_steps = path_to_steps(g, &detour);
                if !detour_steps.iter().any(|s| s.node == target) {
                    continue;
                }

                let mut new_steps: Vec<Step> = walk.steps[..=i].to_vec();
                new_steps.extend(detour_steps);
                if j + 1 < walk.steps.len() {
                    new_steps.extend(walk.steps[(j + 1)..].iter().cloned());
                }
                let new_walk = Walk { steps: new_steps };
                if !is_valid_walk(g, &new_walk) {
                    continue;
                }
                return Some(new_walk);
            }
        }
    }

    None
}

/// Segment-restricted splice: same logic as try_splice but only scans (i, j)
/// pairs within [i_lo, i_hi]. Creates its own BfsCache so it can run on a
/// separate thread without sharing state with other segment tasks.
fn try_splice_in_range(
    g: &Graph,
    walk: &Walk,
    target: usize,
    max_extra_revisits: usize,
    i_lo: usize,
    i_hi: usize,
) -> Option<Walk> {
    if i_hi <= i_lo {
        return None;
    }
    let n_steps = walk.steps.len();
    if n_steps < 2 || i_lo >= n_steps || i_hi >= n_steps {
        return None;
    }

    let original_nodes: HashSet<usize> = walk.steps.iter().map(|s| s.node).collect();
    let mut bfs_cache: BfsCache = HashMap::new();

    let reachable_j: [HashSet<usize>; 2] = [0u8, 1u8].map(|te| {
        let parents_from_t = get_bfs(g, &mut bfs_cache, target, 1 - te);
        ((i_lo + 1)..=i_hi)
            .filter(|&j| {
                let (b_node, b_entry) = entry_state_of(walk, j);
                parents_from_t[halfedge(b_node, 1 - b_entry)].is_some()
            })
            .collect::<HashSet<usize>>()
    });

    if reachable_j[0].is_empty() && reachable_j[1].is_empty() {
        return None;
    }

    let range_len = i_hi - i_lo;
    for gap in 1..=range_len {
        for i in i_lo..(i_hi + 1 - gap) {
            let j = i + gap;
            for target_entry in [0u8, 1u8] {
                if !reachable_j[target_entry as usize].contains(&j) {
                    continue;
                }
                let (a_node, a_exit) = exit_state_after(walk, i);
                let parents_fwd = get_bfs(g, &mut bfs_cache, a_node, a_exit);
                let target_exit_state = halfedge(target, 1 - target_entry);
                if parents_fwd[target_exit_state].is_none() {
                    continue;
                }

                let path_to_target = match reconstruct_path(&parents_fwd, target_exit_state) {
                    Some(p) => p,
                    None => continue,
                };
                let parents_from_t = get_bfs(g, &mut bfs_cache, target, 1 - target_entry);
                let (b_node, b_entry) = entry_state_of(walk, j);
                let target_he = halfedge(b_node, 1 - b_entry);
                let path_from_target = match reconstruct_path(&parents_from_t, target_he) {
                    Some(p) => p,
                    None => continue,
                };

                let mut seen: HashSet<usize> = HashSet::new();
                let mut extras = 0usize;
                for &(n, _, _, _) in path_to_target[1..]
                    .iter()
                    .chain(path_from_target[1..].iter())
                {
                    if original_nodes.contains(&n) || !seen.insert(n) {
                        extras += 1;
                    }
                }
                if extras > max_extra_revisits {
                    continue;
                }

                let mut detour: Vec<(usize, Side, u64, String)> = path_to_target;
                detour.extend(path_from_target.into_iter().skip(1));
                let detour_steps = path_to_steps(g, &detour);
                if !detour_steps.iter().any(|s| s.node == target) {
                    continue;
                }

                let mut new_steps: Vec<Step> = walk.steps[..=i].to_vec();
                new_steps.extend(detour_steps);
                if j + 1 < walk.steps.len() {
                    new_steps.extend(walk.steps[(j + 1)..].iter().cloned());
                }
                let new_walk = Walk { steps: new_steps };
                if !is_valid_walk(g, &new_walk) {
                    continue;
                }
                return Some(new_walk);
            }
        }
    }
    None
}

/// Try to extend a walk past its current end to pick up an uncovered node.
/// This is a simpler variant of try_splice that only appends to the end of
/// the walk. Useful when the walk currently terminates at an internal node
/// because the algorithm decided to stop.
fn try_extend_end(
    g: &Graph,
    walk: &Walk,
    target: usize,
    max_extra_revisits: usize,
    bfs_cache: &mut BfsCache,
) -> Option<Walk> {
    let n_steps = walk.steps.len();
    if n_steps < 1 {
        return None;
    }
    let last = &walk.steps[n_steps - 1];
    let (_, exit) = entry_exit_for_orient(last.orient);
    let parents = get_bfs(g, bfs_cache, last.node, exit);
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
/// into one of the two walks via try_splice or try_extend_end. Splices are
/// only accepted if they strictly increase joint bp coverage.
///
/// The component is partitioned into bubble segments at anchor (cut-vertex)
/// nodes. Within each segment the splice search is restricted to the walk
/// steps between that segment's bounding anchors — reducing the O(L²) scan
/// to O(segment_len²) — and segments are searched in parallel. Uncovered
/// nodes not assigned to any segment fall back to the original sequential
/// full-walk scan.
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
    let mut iters = 0usize;
    let max_iters = members.len() * 4 + 16;

    // Compute anchor structure once for the lifetime of repair.
    // Anchors are the cut-vertex nodes every walk must pass through; they
    // define the segment boundaries for the parallel splice search.
    let member_set: HashSet<usize> = members.iter().copied().collect();
    let anchor_set = find_anchors(g, members, &member_set);
    let ordered: Vec<usize> = if anchor_set.len() >= 2 {
        order_anchors(g, members, &anchor_set, &member_set)
    } else {
        Vec::new()
    };
    let n_segs = ordered.len().saturating_sub(1);

    // Pre-compute the node set for each segment (between consecutive anchors).
    let seg_node_sets: Vec<HashSet<usize>> = (0..n_segs)
        .map(|k| {
            collect_subgraph_nodes(g, &member_set, ordered[k], ordered[k + 1])
                .into_iter()
                .collect()
        })
        .collect();

    // Map each node to the first segment that contains it.
    let node_to_seg: HashMap<usize, usize> = {
        let mut m: HashMap<usize, usize> = HashMap::new();
        for (k, set) in seg_node_sets.iter().enumerate() {
            for &n in set {
                m.entry(n).or_insert(k);
            }
        }
        m
    };

    // Fallback BFS cache for the sequential path (unsegmented nodes, extend_end).
    let mut bfs_cache: BfsCache = HashMap::new();

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

        let cov_before_bp: u64 = covered_before.iter().map(|&n| g.segs[n].length).sum();

        // Partition uncovered nodes: segmented (fast parallel path) vs. the rest.
        let mut by_seg: Vec<Vec<usize>> = vec![Vec::new(); n_segs];
        let mut unsegmented: Vec<usize> = Vec::new();
        for &u in &uncovered {
            match node_to_seg.get(&u) {
                Some(&k) => by_seg[k].push(u),
                None => unsegmented.push(u),
            }
        }

        // Precompute anchor step positions in each walk for this iteration.
        // Must recompute each iteration because prior splices shift step indices.
        let anchor_pos: [HashMap<usize, usize>; 2] = [&w1, &w2].map(|w| {
            ordered
                .iter()
                .filter_map(|&a| w.steps.iter().position(|s| s.node == a).map(|p| (a, p)))
                .collect()
        });
        // For segment k of walk `which`, return (i_lo, i_hi) step bounds or None.
        let seg_range = |which: usize, k: usize| -> Option<(usize, usize)> {
            if k + 1 >= ordered.len() {
                return None;
            }
            let pos = &anchor_pos[which];
            let p = *pos.get(&ordered[k])?;
            let q = *pos.get(&ordered[k + 1])?;
            if p < q { Some((p, q)) } else { None }
        };

        // Precompute step ranges into plain Vecs so the parallel closure captures
        // only Send data (not the seg_range closure which references non-Send locals).
        let w1_ranges: Vec<Option<(usize, usize)>> =
            (0..n_segs).map(|k| seg_range(0, k)).collect();
        let w2_ranges: Vec<Option<(usize, usize)>> =
            (0..n_segs).map(|k| seg_range(1, k)).collect();

        // === Parallel segment splice search ===
        // Each task handles one segment: tries uncovered nodes in that segment
        // (in order) until a splice is found. Returns (which_walk, i_lo, new_walk).
        let candidates: Vec<(usize, usize, Walk)> = (0..n_segs)
            .into_par_iter()
            .filter_map(|k| {
                if by_seg[k].is_empty() {
                    return None;
                }
                let ranges = [w1_ranges[k], w2_ranges[k]];
                for &u in &by_seg[k] {
                    for which in 0..2usize {
                        if let Some((i_lo, i_hi)) = ranges[which] {
                            let walk = if which == 0 { &w1 } else { &w2 };
                            if let Some(new_w) = try_splice_in_range(
                                g, walk, u, max_extra_revisits_per_node, i_lo, i_hi,
                            ) {
                                return Some((which, i_lo, new_w));
                            }
                        }
                    }
                }
                None
            })
            .collect();

        // === Apply candidates ===
        // Separate by which walk, sort descending by i_lo so that when multiple
        // splices land on the same walk the highest-indexed one is applied first
        // (preserving the step-index validity of lower-indexed splices).
        let mut applied = false;

        let apply_cands = |walk: &mut Walk,
                           other: &Walk,
                           cands: &mut Vec<(usize, Walk)>,
                           splices: &mut usize,
                           applied: &mut bool,
                           verbose: bool| {
            cands.sort_by(|a, b| b.0.cmp(&a.0));
            for (_, new_w) in cands.drain(..) {
                let joint: HashSet<usize> =
                    new_w.node_set().union(&other.node_set()).copied().collect();
                let joint_bp: u64 = joint.iter().map(|&n| g.segs[n].length).sum();
                if joint_bp > cov_before_bp {
                    if verbose {
                        eprintln!(
                            "    segmented splice into w (bp {} -> {}, nodes {})",
                            cov_before_bp, joint_bp, joint.len()
                        );
                    }
                    *walk = new_w;
                    *splices += 1;
                    *applied = true;
                }
            }
        };

        let mut w1_cands: Vec<(usize, Walk)> = candidates
            .iter()
            .filter(|(which, _, _)| *which == 0)
            .map(|(_, i_lo, w)| (*i_lo, w.clone()))
            .collect();
        let mut w2_cands: Vec<(usize, Walk)> = candidates
            .iter()
            .filter(|(which, _, _)| *which == 1)
            .map(|(_, i_lo, w)| (*i_lo, w.clone()))
            .collect();

        // Apply w1 splices, then w2 splices. For the w2 coverage check we use
        // the already-updated w1 (monotone: can only help, never hurt).
        apply_cands(&mut w1, &w2, &mut w1_cands, &mut splices, &mut applied, verbose);
        apply_cands(&mut w2, &w1, &mut w2_cands, &mut splices, &mut applied, verbose);

        // === Fallback: sequential full-walk scan for unsegmented nodes ===
        if !applied {
            'splice_search: for u in &unsegmented {
                for which in 0..2usize {
                    let w_ref = if which == 0 { &w1 } else { &w2 };
                    let new_w_opt =
                        try_splice(g, w_ref, *u, max_extra_revisits_per_node, &mut bfs_cache)
                            .or_else(|| {
                                try_extend_end(
                                    g, w_ref, *u, max_extra_revisits_per_node, &mut bfs_cache,
                                )
                            });
                    if let Some(new_w) = new_w_opt {
                        let other = if which == 0 { &w2 } else { &w1 };
                        let joint: HashSet<usize> =
                            new_w.node_set().union(&other.node_set()).copied().collect();
                        let joint_bp: u64 = joint.iter().map(|&n| g.segs[n].length).sum();
                        if joint_bp > cov_before_bp {
                            if verbose {
                                eprintln!(
                                    "    splice {} into w{} (bp {} -> {}, nodes {})",
                                    g.segs[*u].name,
                                    which + 1,
                                    cov_before_bp,
                                    joint_bp,
                                    joint.len(),
                                );
                            }
                            if which == 0 { w1 = new_w; } else { w2 = new_w; }
                            splices += 1;
                            applied = true;
                            break 'splice_search;
                        }
                    }
                }
            }
        }

        // === Fallback: try_extend_end for any still-uncovered segmented nodes ===
        // Handles nodes whose segment anchors weren't found in the walk.
        if !applied {
            'ext_search: for &u in &uncovered {
                for which in 0..2usize {
                    let w_ref = if which == 0 { &w1 } else { &w2 };
                    if let Some(new_w) =
                        try_extend_end(g, w_ref, u, max_extra_revisits_per_node, &mut bfs_cache)
                    {
                        let other = if which == 0 { &w2 } else { &w1 };
                        let joint: HashSet<usize> =
                            new_w.node_set().union(&other.node_set()).copied().collect();
                        let joint_bp: u64 = joint.iter().map(|&n| g.segs[n].length).sum();
                        if joint_bp > cov_before_bp {
                            if which == 0 { w1 = new_w; } else { w2 = new_w; }
                            splices += 1;
                            applied = true;
                            break 'ext_search;
                        }
                    }
                }
            }
        }

        if !applied {
            if verbose {
                eprintln!(
                    "    coverage repair stalled with {} node(s) still uncovered",
                    uncovered.len()
                );
            }
            break;
        }
    }

    let final_uncovered = members
        .iter()
        .filter(|u| !w1.node_set().contains(u) && !w2.node_set().contains(u))
        .count();
    (w1, w2, splices, final_uncovered)
}

// ------------ FASTA helpers ------------

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
            "usage: {} <input.gfa> [--algorithm NAME] [--endpoints tip|any] [--candidates K] [--seed S] [--sample NAME] [--skip-sequences] [--length-tolerance R] [--repair-max-revisits N] [--no-repair] [--threads N] [--verbose]",
            args.get(0).map(|s| s.as_str()).unwrap_or("gfa_haps")
        );
        process::exit(2);
    }

    let input_path = PathBuf::from(&args[1]);
    let mut algo_name = String::from("diverse-pair");
    let mut mode = EndpointMode::Tip;
    let mut k: usize = 64;
    let mut seed: u64 = 0xC0FFEE;
    let mut sample = String::from("sample");
    let mut skip_sequences = false;
    let mut verbose = false;
    let mut repair_max_revisits: usize = 4;
    let mut no_repair = false;
    let mut length_tolerance: f64 = 0.10;
    let mut threads: usize = 0; // 0 = let rayon choose (all available cores)

    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--algorithm" => {
                i += 1;
                if i >= args.len() {
                    die("--algorithm needs a value");
                }
                algo_name = args[i].clone();
            }
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
            "--length-tolerance" => {
                i += 1;
                if i >= args.len() {
                    die("--length-tolerance needs a value");
                }
                length_tolerance = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--length-tolerance must be a number in [0, 1]"));
                if length_tolerance < 0.0 || length_tolerance > 1.0 {
                    die("--length-tolerance must be between 0.0 and 1.0");
                }
            }
            "--threads" => {
                i += 1;
                if i >= args.len() {
                    die("--threads needs a value");
                }
                threads = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--threads must be a positive integer"));
                if threads == 0 {
                    die("--threads must be >= 1");
                }
            }
            other => die(&format!("unknown argument: {}", other)),
        }
        i += 1;
    }

    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap_or_else(|e| die(&format!("failed to build thread pool: {}", e)));
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

    let algorithm: Box<dyn WalkPairAlgorithm> = match algo_name.as_str() {
        "diverse-pair" => Box::new(DiversePair { mode, k, seed, max_length_ratio: length_tolerance }),
        "divide-and-conquer" => Box::new(DivideAndConquer { mode, k, seed, max_length_ratio: length_tolerance }),
        other => die(&format!("unknown algorithm: {}", other)),
    };

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

        match algorithm.build_pair(g, m, ci, verbose) {
            None => {
                eprintln!("warning: no walks generated for {}", component_label);
            }
            Some((w1, w2)) => {
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
                let sc: PairScore = pair_score(g, &w1, &w2);
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
