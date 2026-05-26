// divide_and_conquer: partition a component at mandatory anchor nodes,
// solve each bubble subgraph independently, then chain the sub-solutions.
//
// Stages:
//   1. Find anchor nodes — nodes whose removal disconnects the component's
//      two reference tips from each other. Every full walk must traverse
//      each anchor.
//   2. Order anchors linearly via BFS from a tip anchor.
//   3. For each consecutive anchor pair, collect the subgraph via two-BFS
//      intersection: nodes reachable from a_start (not past a_end) ∩ nodes
//      reachable from a_end (not past a_start).
//   4. Solve each subgraph with best_pair_for_component (node_mask confined).
//      Trim each walk to [a_start .. a_end]. Walks in the wrong direction
//      are detected and reversed with correct orientation/overlap fixup. If
//      a walk doesn't reach a_end at all (e.g. pair-scoring picked a short
//      walk to maximise symdiff), replace it with a BFS-derived anchor walk.
//   5. Chain sub-solutions with greedy DP: at each junction try same/flipped
//      phase assignment, keep the better pair_score.
//
// Falls back to DiversePair on the whole component whenever decomposition
// cannot produce a valid chain (too few anchors, failed trim, incompatible
// junction orientation, etc.).

use gfaphaser::{Graph, Step, Walk};
use std::collections::{HashSet, VecDeque};
use super::WalkPairAlgorithm;
use super::diverse_pair::{best_pair_for_component, pair_better, pair_score, EndpointMode};

// ------------ BFS helpers ------------

/// Undirected BFS from `start` within `member_set`. Visits `barrier` but
/// does not expand through it, so nodes "beyond" the barrier are excluded.
fn bfs_reachable_set(
    g: &Graph,
    start: usize,
    barrier: usize,
    member_set: &HashSet<usize>,
) -> HashSet<usize> {
    let mut visited: HashSet<usize> = HashSet::new();
    let mut q: VecDeque<usize> = VecDeque::new();
    visited.insert(start);
    q.push_back(start);
    while let Some(cur) = q.pop_front() {
        if cur == barrier {
            continue;
        }
        for side in [0u8, 1u8] {
            for e in g.neighbors(cur, side) {
                if member_set.contains(&e.to) && visited.insert(e.to) {
                    q.push_back(e.to);
                }
            }
        }
    }
    visited
}

/// Is `to` reachable from `from` within `member_set` when `skip` is removed?
fn bfs_reachable(
    g: &Graph,
    from: usize,
    to: usize,
    member_set: &HashSet<usize>,
    skip: usize,
) -> bool {
    if from == skip || to == skip {
        return false;
    }
    let mut visited: HashSet<usize> = HashSet::new();
    let mut q: VecDeque<usize> = VecDeque::new();
    visited.insert(from);
    q.push_back(from);
    while let Some(cur) = q.pop_front() {
        if cur == to {
            return true;
        }
        for side in [0u8, 1u8] {
            for e in g.neighbors(cur, side) {
                if e.to != skip && member_set.contains(&e.to) && visited.insert(e.to) {
                    q.push_back(e.to);
                }
            }
        }
    }
    false
}

// ------------ Stage 1: find anchors ------------

/// Return the set of anchor nodes: every node whose removal disconnects
/// the component's reference tips (tip_a, tip_b). Both tips are included.
/// Returns an empty set when there are fewer than 2 tips (no useful
/// decomposition possible).
fn find_anchors(g: &Graph, members: &[usize], member_set: &HashSet<usize>) -> HashSet<usize> {
    let tips: Vec<usize> = members
        .iter()
        .copied()
        .filter(|&u| g.is_tip_on(u, 0) || g.is_tip_on(u, 1))
        .collect();
    if tips.len() < 2 {
        return HashSet::new();
    }
    let tip_a = tips[0];
    let tip_b = tips[tips.len() - 1];
    let mut anchors = HashSet::new();
    anchors.insert(tip_a);
    anchors.insert(tip_b);
    for &n in members {
        if n == tip_a || n == tip_b {
            continue;
        }
        if !bfs_reachable(g, tip_a, tip_b, member_set, n) {
            anchors.insert(n);
        }
    }
    anchors
}

// ------------ Stage 2: order anchors ------------

/// BFS from a tip anchor; record anchors in the order they are first visited.
fn order_anchors(
    g: &Graph,
    members: &[usize],
    anchors: &HashSet<usize>,
    member_set: &HashSet<usize>,
) -> Vec<usize> {
    let start = match members
        .iter()
        .copied()
        .find(|&u| anchors.contains(&u) && (g.is_tip_on(u, 0) || g.is_tip_on(u, 1)))
    {
        Some(s) => s,
        None => return Vec::new(),
    };
    let mut visited: HashSet<usize> = HashSet::new();
    let mut q: VecDeque<usize> = VecDeque::new();
    let mut ordered: Vec<usize> = Vec::new();
    visited.insert(start);
    q.push_back(start);
    ordered.push(start); // start is always an anchor (it's a tip)
    while let Some(cur) = q.pop_front() {
        for side in [0u8, 1u8] {
            for e in g.neighbors(cur, side) {
                if member_set.contains(&e.to) && visited.insert(e.to) {
                    if anchors.contains(&e.to) {
                        ordered.push(e.to);
                    }
                    q.push_back(e.to);
                }
            }
        }
    }
    ordered
}

// ------------ Stage 3: collect subgraph nodes ------------

/// Nodes that lie between a_start and a_end: intersection of the set
/// reachable from a_start (without passing through a_end) and the set
/// reachable from a_end (without passing through a_start).
fn collect_subgraph_nodes(
    g: &Graph,
    member_set: &HashSet<usize>,
    a_start: usize,
    a_end: usize,
) -> Vec<usize> {
    let fwd = bfs_reachable_set(g, a_start, a_end, member_set);
    let bwd = bfs_reachable_set(g, a_end, a_start, member_set);
    fwd.intersection(&bwd).copied().collect()
}

// ------------ Walk trimming and reversal ------------

/// Reverse a walk, flipping each step's orientation and adjusting per-step
/// overlap_in/cigar_in values for the reversed direction.
///
/// For a forward walk [s0, s1, ..., s_{n-1}]:
///   reversed = [s'_{n-1}, ..., s'_0]  (orient flipped)
///   rev[k].overlap_in = original s_{n-k}.overlap_in  (k >= 1), else 0
///
/// Both directions of an L-line share the same CIGAR overlap, so this
/// assignment is correct.
fn reverse_walk(walk: Walk) -> Walk {
    if walk.steps.is_empty() {
        return walk;
    }
    let mut steps: Vec<_> = walk
        .steps
        .into_iter()
        .rev()
        .map(|mut s| {
            s.orient = if s.orient == '+' { '-' } else { '+' };
            s
        })
        .collect();
    // Rotate overlaps left by one position, inserting 0 at index 0.
    // After the Vec reversal, steps[k] holds original step[n-1-k].
    // We want steps[k].overlap_in = original step[n-k].overlap_in = initial steps[k-1].overlap_in.
    let mut prev_ov = 0u64;
    let mut prev_cig = String::new();
    for step in &mut steps {
        let this_ov = step.overlap_in;
        let this_cig = step.cigar_in.clone();
        step.overlap_in = prev_ov;
        step.cigar_in = prev_cig;
        prev_ov = this_ov;
        prev_cig = this_cig;
    }
    Walk { steps }
}

/// Trim a walk to the range [first a_start .. first a_end after a_start].
/// If the walk runs in the reverse direction (a_end before a_start), the
/// relevant segment is extracted and reversed. Returns None only if neither
/// endpoint appears in the walk.
fn trim_walk(walk: Walk, a_start: usize, a_end: usize) -> Option<Walk> {
    // Forward direction: a_start → a_end
    let si = walk.steps.iter().position(|s| s.node == a_start);
    if let Some(si) = si {
        if let Some(rel_ei) = walk.steps[si..].iter().position(|s| s.node == a_end) {
            return Some(Walk {
                steps: walk.steps[si..=si + rel_ei].to_vec(),
            });
        }
    }
    // Reversed direction: a_end → a_start, then reverse the trimmed walk
    let si_r = walk.steps.iter().position(|s| s.node == a_end)?;
    let rel_ei_r = walk.steps[si_r..].iter().position(|s| s.node == a_start)?;
    let trimmed_rev = Walk {
        steps: walk.steps[si_r..=si_r + rel_ei_r].to_vec(),
    };
    Some(reverse_walk(trimmed_rev))
}

/// Bidirected-state BFS: find any valid walk from `start_node` to
/// `target_node` confined to `sub_set`. Used as a fallback when the walk
/// returned by best_pair_for_component does not reach the endpoint anchor.
fn bfs_anchor_walk(
    g: &Graph,
    start_node: usize,
    target_node: usize,
    sub_set: &HashSet<usize>,
) -> Option<Walk> {
    let n = g.n();
    // parent[halfedge(v, exit)] = (parent_halfedge, overlap, cigar) or None
    let mut parent: Vec<Option<(usize, u64, String)>> = vec![None; n * 2];
    let mut q: VecDeque<usize> = VecDeque::new();

    // Seed both exit sides of start_node
    for start_exit in [1u8, 0u8] {
        let he = start_node * 2 + start_exit as usize;
        if parent[he].is_none() {
            parent[he] = Some((usize::MAX, 0, String::new()));
            q.push_back(he);
        }
    }

    let mut found_he: Option<usize> = None;
    'bfs: while let Some(h) = q.pop_front() {
        let node = h / 2;
        if node == target_node {
            found_he = Some(h);
            break 'bfs;
        }
        let exit_side = (h % 2) as u8;
        for e in g.neighbors(node, exit_side) {
            if !sub_set.contains(&e.to) {
                continue;
            }
            let next_exit = 1 - e.to_side;
            let nh = e.to * 2 + next_exit as usize;
            if parent[nh].is_none() {
                parent[nh] = Some((h, e.overlap, e.cigar.clone()));
                q.push_back(nh);
            }
        }
    }

    let target_he = found_he?;

    // Reconstruct the path of halfedges from start to target
    let mut he_path: Vec<(usize, u64, String)> = Vec::new();
    let mut cur = target_he;
    loop {
        let entry = parent[cur].as_ref().unwrap();
        he_path.push((cur, entry.1, entry.2.clone()));
        if entry.0 == usize::MAX {
            break;
        }
        cur = entry.0;
    }
    he_path.reverse();

    // Convert halfedge path to Walk steps
    let steps: Vec<Step> = he_path
        .iter()
        .enumerate()
        .map(|(i, (he, ov, cig))| {
            let node = he / 2;
            let exit_side = (he % 2) as u8;
            let entry_side = 1 - exit_side;
            let orient = if entry_side == 0 { '+' } else { '-' };
            Step {
                node,
                orient,
                overlap_in: if i == 0 { 0 } else { *ov },
                cigar_in: if i == 0 { String::new() } else { cig.clone() },
            }
        })
        .collect();

    if steps.is_empty() {
        return None;
    }
    Some(Walk { steps })
}

// ------------ Stage 5 helpers ------------

/// Concatenate two walks, dropping the duplicate anchor step at the junction.
/// `left.steps.last()` and `right.steps[0]` are the same anchor node; the
/// duplicate is dropped from the right walk's prefix.
fn concat_walks(left: &Walk, right: &Walk) -> Walk {
    let mut steps = left.steps.clone();
    if right.steps.len() > 1 {
        steps.extend_from_slice(&right.steps[1..]);
    }
    Walk { steps }
}

// ------------ Algorithm implementation ------------

pub struct DivideAndConquer {
    pub mode: EndpointMode,
    pub k: usize,
    pub seed: u64,
    pub max_length_ratio: f64,
}

impl WalkPairAlgorithm for DivideAndConquer {
    fn build_pair(
        &self,
        g: &Graph,
        members: &[usize],
        comp_idx: usize,
        verbose: bool,
    ) -> Option<(Walk, Walk)> {
        let comp_seed = self
            .seed
            .wrapping_add((comp_idx as u64).wrapping_mul(0x9E3779B97F4A7C15));

        macro_rules! fallback {
            () => {
                best_pair_for_component(
                    g, members, self.mode, self.k, comp_seed, None, verbose, self.max_length_ratio,
                )
            };
        }

        let member_set: HashSet<usize> = members.iter().copied().collect();

        // Stage 1 — anchors
        let anchor_set = find_anchors(g, members, &member_set);
        if anchor_set.len() < 2 {
            if verbose {
                eprintln!("  divide-and-conquer: no anchors, falling back");
            }
            return fallback!();
        }

        // Stage 2 — linear ordering
        let ordered = order_anchors(g, members, &anchor_set, &member_set);
        if ordered.len() <= 2 {
            if verbose {
                eprintln!("  divide-and-conquer: no interior anchors, falling back");
            }
            return fallback!();
        }

        if verbose {
            eprintln!(
                "  divide-and-conquer: {} anchors, {} subgraphs for {} nodes",
                ordered.len(),
                ordered.len() - 1,
                members.len()
            );
        }

        // Stages 3 & 4 — collect and solve each subgraph
        let mut solutions: Vec<(Walk, Walk)> = Vec::new();
        let mut subgraph_seed = comp_seed;

        for i in 0..(ordered.len() - 1) {
            let a_start = ordered[i];
            let a_end = ordered[i + 1];

            let sub_nodes = collect_subgraph_nodes(g, &member_set, a_start, a_end);
            let sub_set: HashSet<usize> = sub_nodes.iter().copied().collect();

            subgraph_seed = subgraph_seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);

            let (raw_w1, raw_w2) = match best_pair_for_component(
                g,
                &sub_nodes,
                EndpointMode::Any,
                self.k,
                subgraph_seed,
                Some(&sub_set),
                false, // suppress per-subgraph verbose pool dumps
                1.0,   // no size filter for internal subgraph calls
            ) {
                Some(pair) => pair,
                None => {
                    if verbose {
                        eprintln!(
                            "  divide-and-conquer: no solution for subgraph {}->{}, falling back",
                            a_start, a_end
                        );
                    }
                    return fallback!();
                }
            };

            // Trim each walk to [a_start..a_end].  If trim fails (walk doesn't
            // reach a_end), substitute the BFS shortest-path anchor walk.
            let w1 = trim_walk(raw_w1, a_start, a_end)
                .or_else(|| bfs_anchor_walk(g, a_start, a_end, &sub_set));
            let w2 = trim_walk(raw_w2, a_start, a_end)
                .or_else(|| bfs_anchor_walk(g, a_start, a_end, &sub_set));

            match (w1, w2) {
                (Some(w1), Some(w2)) => solutions.push((w1, w2)),
                _ => {
                    if verbose {
                        eprintln!(
                            "  divide-and-conquer: no anchor walk found for {}->{}, falling back",
                            a_start, a_end
                        );
                    }
                    return fallback!();
                }
            }
        }

        if solutions.is_empty() {
            return fallback!();
        }

        // Stage 5 — greedy DP chain
        let (mut h1, mut h2) = (solutions[0].0.clone(), solutions[0].1.clone());

        for sol in &solutions[1..] {
            let (r1, r2) = (&sol.0, &sol.1);

            let l1_orient = h1.steps.last().unwrap().orient;
            let l2_orient = h2.steps.last().unwrap().orient;
            let r1_orient = r1.steps.first().unwrap().orient;
            let r2_orient = r2.steps.first().unwrap().orient;

            // "same" phase: h1 continues with r1, h2 continues with r2
            let same_ok = l1_orient == r1_orient && l2_orient == r2_orient;
            // "flip" phase: h1 continues with r2, h2 continues with r1
            let flip_ok = l1_orient == r2_orient && l2_orient == r1_orient;

            if !same_ok && !flip_ok {
                if verbose {
                    eprintln!(
                        "  divide-and-conquer: junction orientation incompatible, falling back"
                    );
                }
                return fallback!();
            }

            let mut candidates: Vec<(Walk, Walk)> = Vec::with_capacity(2);
            if same_ok {
                candidates.push((concat_walks(&h1, r1), concat_walks(&h2, r2)));
            }
            if flip_ok {
                candidates.push((concat_walks(&h1, r2), concat_walks(&h2, r1)));
            }

            let mut best_idx = 0;
            let mut best_sc = pair_score(g, &candidates[0].0, &candidates[0].1);
            for (idx, (c1, c2)) in candidates.iter().enumerate().skip(1) {
                let sc = pair_score(g, c1, c2);
                if pair_better(sc, best_sc) {
                    best_idx = idx;
                    best_sc = sc;
                }
            }

            let (new_h1, new_h2) = candidates.remove(best_idx);
            h1 = new_h1;
            h2 = new_h2;
        }

        Some((h1, h2))
    }
}
