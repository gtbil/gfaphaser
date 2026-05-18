// diverse_pair: the "diverse-pair" walk-pair algorithm.
//
// Generates K randomized candidate walks, augments the pool with
// "complementary" walks biased to diverge from the best single walk
// (haplotype-divergence guide), then selects the best pair by a
// lexicographic objective:
//   1. maximize union coverage (bp covered by either walk)
//   2. maximize symmetric difference (bp covered by exactly one walk)
//   3. minimize total revisits

use gfaphaser::{Edge, Graph, Side, Step, Walk};
use std::collections::HashSet;
use super::WalkPairAlgorithm;

// ------------ Simple deterministic PRNG (xorshift64*) ------------
//
// Avoids pulling in the `rand` crate so we stay dependency-free.

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
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

// ------------ Endpoint mode ------------

#[derive(Clone, Copy)]
pub enum EndpointMode {
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
    if exit_side == 1 { '+' } else { '-' }
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
    node_mask: Option<&HashSet<usize>>,
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
        let masked_buf: Option<Vec<Edge>> = node_mask.map(|m| {
            g.neighbors(cur, cur_exit)
                .iter()
                .filter(|e| m.contains(&e.to))
                .cloned()
                .collect()
        });
        let edges: &[Edge] = match &masked_buf {
            Some(v) => v.as_slice(),
            None => g.neighbors(cur, cur_exit),
        };
        if edges.is_empty() {
            break;
        }

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
    node_mask: Option<&HashSet<usize>>,
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

        let walk = generate_walk(g, sn, se, &mut rng, budget, None, 1.0, node_mask);
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

// ------------ Pair scoring ------------

#[derive(Clone, Copy, Debug)]
pub struct PairScore {
    pub union_bp: u64,
    pub symdiff_bp: u64,
    pub revisits_total: i64,
}

pub fn pair_score(g: &Graph, w1: &Walk, w2: &Walk) -> PairScore {
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

pub fn pair_better(a: PairScore, b: PairScore) -> bool {
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
pub fn best_pair_for_component(
    g: &Graph,
    members: &[usize],
    mode: EndpointMode,
    k: usize,
    base_seed: u64,
    node_mask: Option<&HashSet<usize>>,
    verbose: bool,
) -> Option<(Walk, Walk)> {
    let mut pool = generate_candidates(g, members, mode, k, base_seed, node_mask);
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
        let w = generate_walk(g, sn, se, &mut rng, budget, Some(&leader_nodes), gw, node_mask);
        pool.push(w);
    }

    // Evaluate all pairs (ordered as the canonical (smaller_idx, larger_idx)).
    if verbose {
        eprintln!("  pool size: {}", pool.len());
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

    best.map(|(i, j, _sc)| (pool[i].clone(), pool[j].clone()))
}

// ------------ Algorithm implementation ------------

/// Randomized candidate generation + lexicographic best-pair selection
/// (max union bp → max symdiff bp → min revisits), with complementary walks
/// biased to diverge through bubbles. Name: "diverse-pair".
pub struct DiversePair {
    pub mode: EndpointMode,
    pub k: usize,
    pub seed: u64,
}

impl WalkPairAlgorithm for DiversePair {
    fn build_pair(&self, g: &Graph, members: &[usize], comp_idx: usize, verbose: bool) -> Option<(Walk, Walk)> {
        let comp_seed = self.seed
            .wrapping_add((comp_idx as u64).wrapping_mul(0x9E3779B97F4A7C15));
        best_pair_for_component(g, members, self.mode, self.k, comp_seed, None, verbose)
    }
}
