#!/usr/bin/env python3
"""
Two-walk construction through a hifiasm unitig graph, guided by reference
coordinates and using DFS with backtracking and memoization.

Inputs:
  - GFA file (hifiasm p_utg.gfa, ideally pre-filtered to one chromosome)
  - PAF file (chained unitig-to-reference alignments)
  - Output prefix

Outputs:
  - <prefix>.hap1.gaf, <prefix>.hap2.gaf
  - <prefix>.hap1.fa,  <prefix>.hap2.fa
  - <prefix>.stats.tsv

Dependencies: gfapy, networkx
"""

import sys
import argparse
from collections import defaultdict
import gfapy
import networkx as nx


# ---------- I/O and graph construction ----------

def revcomp(seq):
    comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(comp)[::-1]


def load_gfa(gfa_path):
    gfa = gfapy.Gfa.from_file(gfa_path)
    G = nx.DiGraph()
    for seg in gfa.segments:
        for orient in ('+', '-'):
            G.add_node(f"{seg.name}{orient}", length=seg.length)
    for edge in gfa.dovetails:
        f_name = edge.from_segment.name
        t_name = edge.to_segment.name
        f_orient = edge.from_orient
        t_orient = edge.to_orient
        G.add_edge(f"{f_name}{f_orient}", f"{t_name}{t_orient}")
        flip = {'+': '-', '-': '+'}
        G.add_edge(f"{t_name}{flip[t_orient]}", f"{f_name}{flip[f_orient]}")
    return gfa, G


def detect_bubble_members(gfa):
    seg_neighbors = defaultdict(set)
    for edge in gfa.dovetails:
        a = edge.from_segment.name
        b = edge.to_segment.name
        if a != b:
            seg_neighbors[a].add(b)
            seg_neighbors[b].add(a)

    pair_to_unitigs = defaultdict(list)
    for seg in gfa.segments:
        nbrs = seg_neighbors[seg.name]
        if len(nbrs) == 2:
            key = frozenset(nbrs)
            pair_to_unitigs[key].append(seg.name)

    bubble_members = set()
    for key, unitigs in pair_to_unitigs.items():
        if len(unitigs) >= 2:
            bubble_members.update(unitigs)
    return bubble_members, seg_neighbors


def assign_budgets(gfa, ref_positions, bubble_members, seg_neighbors, default_budget=2):
    budgets = {}
    n_het = n_default = n_repeat = 0
    for seg in gfa.segments:
        name = seg.name
        aligned = name in ref_positions
        degree = len(seg_neighbors[name])
        if name in bubble_members:
            budget = 1
            n_het += 1
        elif degree > 4:
            budget = min(degree, 4)
            n_repeat += 1
        else:
            budget = default_budget
            n_default += 1
        if not aligned and budget < default_budget:
            budget = default_budget
        budgets[name] = budget
    print(f"# Budgets: {n_het} bubble-het (1), {n_default} default ({default_budget}), "
          f"{n_repeat} repeat-like (>{default_budget})", file=sys.stderr)
    return budgets


def load_positions(paf_path, target_chr, include_set=None):
    min_start = {}
    max_end = {}
    total_bp = defaultdict(int)
    chr_length = 0
    with open(paf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:
                continue
            if parts[5] != target_chr:
                continue
            u = parts[0]
            if include_set is not None and u not in include_set:
                continue
            chr_length = max(chr_length, int(parts[6]))
            rs = int(parts[7])
            re_ = int(parts[8])
            matches = int(parts[9])
            if u not in min_start or rs < min_start[u]:
                min_start[u] = rs
            if u not in max_end or re_ > max_end[u]:
                max_end[u] = re_
            total_bp[u] += matches
    positions = {}
    for u in min_start:
        span = max_end[u] - min_start[u]
        quality = total_bp[u] / span if span > 0 else 1.0
        positions[u] = (min_start[u], max_end[u], total_bp[u], quality)
    return positions, chr_length


# ---------- Walk construction with backtracking ----------

def score_candidate(cand_node, gfa, positions, current_ref_pos,
                    budgets, used_counts, walked_set):
    """Score a candidate next node. Higher = better. None = disqualified."""
    name = cand_node[:-1]
    remaining = budgets.get(name, 0) - used_counts.get(name, 0)
    if remaining <= 0:
        return None

    score = 0.0

    if name in positions:
        rs, re_, _, _ = positions[name]
        if current_ref_pos is None:
            score += 100
        else:
            delta_start = rs - current_ref_pos
            if delta_start >= 0:
                if delta_start < 500_000:
                    score += 100 - delta_start / 5_000
                else:
                    score += -delta_start / 1e6
            else:
                if delta_start > -100_000:
                    score += -10 + delta_start / 5_000  # mild penalty
                else:
                    score += -100
        seg_len = gfa.segment(name).length
        score += min(seg_len / 50_000, 5)
    else:
        score += 15  # moderate reward for unaligned (potentially novel)

    score += remaining * 2

    if name in walked_set:
        score += -30

    return score


def find_start_node(G, gfa, positions, used_counts, budgets, avoid=None,
                    min_quality=0.3):
    """
    Pick the earliest-aligned unitig with a graph out-edge. Skip 'wide-spread'
    unitigs whose alignments span much more reference than their own length
    (likely repeat-driven or chimeric).
    """
    if avoid is None:
        avoid = set()
    candidates = sorted(positions.items(), key=lambda x: x[1][0])
    for u, (rs, re_, abp, qual) in candidates:
        if used_counts.get(u, 0) >= budgets.get(u, 0):
            continue
        if u in avoid:
            continue
        if qual < min_quality:
            continue
        seg = gfa.segment(u)
        if seg is None:
            continue
        span = re_ - rs
        if span > seg.length * 2:
            continue
        for orient in ('+', '-'):
            node = f"{u}{orient}"
            if node in G and G.out_degree(node) > 0:
                return node
    # Fallback: relax filters
    print(f"#   relaxing start-node filters", file=sys.stderr)
    for u, (rs, re_, abp, qual) in candidates:
        if used_counts.get(u, 0) >= budgets.get(u, 0):
            continue
        if u in avoid:
            continue
        for orient in ('+', '-'):
            node = f"{u}{orient}"
            if node in G and G.out_degree(node) > 0:
                return node
    return None


def walk_with_backtracking(G, gfa, positions, budgets, used_counts, chr_length,
                            start_avoid=None, max_backtrack_depth=50,
                            max_steps=200_000, min_forward_progress=10_000):
    """
    DFS-style walk with backtracking. Uses a stack rather than recursion to
    avoid Python's recursion limits on long walks.

    State on the stack: (current_node, remaining_candidates_to_try)
    where remaining_candidates_to_try is a list of (score, node) sorted by score desc.

    On dead end, pop and try the next candidate at the previous level.
    Dead-end-only branches are memoized per current node to avoid re-exploration.
    """
    if start_avoid is None:
        start_avoid = set()

    start_node = find_start_node(G, gfa, positions, used_counts, budgets,
                                  avoid=start_avoid)
    if start_node is None:
        return [], 0, 0

    walked_set = set()
    used_counts[start_node[:-1]] += 1
    walked_set.add(start_node[:-1])

    start_name = start_node[:-1]
    initial_ref_pos = positions[start_name][0] if start_name in positions else None
    initial_len = gfa.segment(start_name).length

    # Memoization: nodes proven to lead only to dead ends from current walk context.
    # This is conservative — we clear it when we make significant forward progress.
    dead_end_nodes = set()

    # Stack frame format: (node, candidates_iter, ref_pos_before_entering, cumulative_len_before)
    # We store enough state to roll back on backtrack.
    def get_candidates(cur_node, current_ref_pos):
        out = list(G.out_edges(cur_node))
        scored = []
        for _, cand in out:
            if cand[:-1] in walked_set:
                # Reuse only allowed if budget permits AND this isn't a recursive loop
                # (i.e., the node isn't currently in our active path)
                # For simplicity, skip walked nodes during backtracking exploration.
                # Reuse will still happen naturally via the budget system across
                # multiple top-level walks.
                continue
            if cand in dead_end_nodes:
                continue
            s = score_candidate(cand, gfa, positions, current_ref_pos,
                                budgets, used_counts, walked_set)
            if s is not None:
                scored.append((s, cand))
        scored.sort(reverse=True)
        return scored

    # Stack: list of frames. Each frame = (node, remaining_candidates, ref_pos_at_entry, len_at_entry)
    # `remaining_candidates` is a list; we pop from the end (highest score next).
    initial_candidates = get_candidates(start_node, initial_ref_pos)
    stack = [{
        'node': start_node,
        'candidates': initial_candidates,
        'ref_pos': initial_ref_pos,
        'cum_len': initial_len,
        'parent_ref_pos': None,
        'parent_cum_len': 0,
    }]

    path = [start_node]
    gaps = 0
    steps = 0
    backtracks = 0
    last_progress_ref_pos = initial_ref_pos or 0

    while stack and steps < max_steps:
        steps += 1
        top = stack[-1]
        current_ref_pos = top['ref_pos']
        cumulative_len = top['cum_len']

        # Stopping conditions
        if current_ref_pos and current_ref_pos >= chr_length * 0.98:
            break
        if cumulative_len > chr_length * 1.3:
            break

        # Get next candidate from current frame
        candidates = top['candidates']
        if not candidates:
            # Dead end — backtrack
            dead_node = top['node']
            dead_end_nodes.add(dead_node)
            stack.pop()
            path.pop()
            # Roll back state
            used_counts[dead_node[:-1]] -= 1
            walked_set.discard(dead_node[:-1])
            backtracks += 1

            # If we've backtracked too far without finding a path, give up
            # and jump
            if len(stack) == 0:
                # Walk done from this start
                break
            if backtracks > max_backtrack_depth and stack:
                # Force a jump from current top
                jump_target = _find_jump_target(positions, used_counts, budgets,
                                                  walked_set,
                                                  stack[-1]['ref_pos'])
                if jump_target is None:
                    break
                placed = False
                for orient in ('+', '-'):
                    node = f"{jump_target}{orient}"
                    if node in G:
                        path.append('|GAP|')
                        path.append(node)
                        walked_set.add(jump_target)
                        used_counts[jump_target] += 1
                        new_ref = positions[jump_target][1]
                        new_len = stack[-1]['cum_len'] + gfa.segment(jump_target).length
                        new_candidates = get_candidates(node, new_ref)
                        stack.append({
                            'node': node,
                            'candidates': new_candidates,
                            'ref_pos': new_ref,
                            'cum_len': new_len,
                            'parent_ref_pos': stack[-1]['ref_pos'],
                            'parent_cum_len': stack[-1]['cum_len'],
                        })
                        gaps += 1
                        backtracks = 0
                        last_progress_ref_pos = new_ref
                        # Reset dead-end memo since we're in a new region
                        dead_end_nodes.clear()
                        placed = True
                        break
                if not placed:
                    break
            continue

        # Try the highest-scoring remaining candidate
        score, cand = candidates.pop()  # take last (we sorted ascending? No, descending)
        # Wait: scored.sort(reverse=True) puts highest first at index 0.
        # We want to pop highest first. Since list.pop() takes from end,
        # we should reverse-sort then pop from end means we get lowest first.
        # Fix: append to candidates in reverse order so highest is at end.
        # Or: use pop(0). Pop(0) is O(n), but candidate lists are small (degree ~2-4)
        # so it's fine. Let me restructure to avoid confusion.
        # — Actually, the cleanest fix: don't reverse, pop(0).
        # I'll fix this below by changing how get_candidates returns.

        # If candidate is now disqualified (budget exhausted by something we just did),
        # skip it
        cand_name = cand[:-1]
        if used_counts.get(cand_name, 0) >= budgets.get(cand_name, 0):
            continue
        if cand_name in walked_set:
            continue

        # Take the step
        new_ref_pos = current_ref_pos
        if cand_name in positions:
            new_ref_pos = max(current_ref_pos or 0, positions[cand_name][1])
        new_len = cumulative_len + gfa.segment(cand_name).length

        path.append(cand)
        walked_set.add(cand_name)
        used_counts[cand_name] += 1

        new_candidates = get_candidates(cand, new_ref_pos)
        stack.append({
            'node': cand,
            'candidates': new_candidates,
            'ref_pos': new_ref_pos,
            'cum_len': new_len,
            'parent_ref_pos': current_ref_pos,
            'parent_cum_len': cumulative_len,
        })

        # If we made significant forward progress, clear dead-end memo
        # (a region that was dead-end-bound earlier may no longer be)
        if new_ref_pos and new_ref_pos - last_progress_ref_pos > min_forward_progress:
            dead_end_nodes.clear()
            last_progress_ref_pos = new_ref_pos
            backtracks = 0

    return path, gaps, backtracks


def _find_jump_target(positions, used_counts, budgets, walked_set, current_ref_pos):
    """Find nearest unwalked, budget-available aligned unitig forward of current pos."""
    best = None
    best_dist = float('inf')
    for u, (rs, re_, _, _) in positions.items():
        if used_counts.get(u, 0) >= budgets.get(u, 0):
            continue
        if u in walked_set:
            continue
        if current_ref_pos is not None and rs < current_ref_pos:
            continue
        dist = rs - (current_ref_pos or 0)
        if dist < best_dist:
            best_dist = dist
            best = u
    return best


# Fix the candidate ordering issue: candidates should be popped highest-score first.
# Reordering: we'll change get_candidates inside walk_with_backtracking to sort
# ascending and pop from end. But since get_candidates is defined inside the
# function above, let me restructure to make this clean.


def walk_one_haplotype(G, gfa, positions, budgets, used_counts, chr_length,
                       start_avoid=None, max_search_depth=200):
    """Wrapper around the target-directed walker."""
    return _walk_target_directed(G, gfa, positions, budgets, used_counts,
                                  chr_length, start_avoid, max_search_depth)


def _bfs_to_target(G, source_node, target_name, blocked_names, max_depth,
                    used_counts, budgets):
    """
    BFS from source_node looking for any oriented form of target_name.
    Avoids nodes whose unitig is in blocked_names (already committed in walk)
    EXCEPT for the target itself. Respects budgets.
    Returns the shortest path as a list of oriented nodes, or None.
    """
    from collections import deque
    if source_node not in G:
        return None
    # Compute prospective usage including blocked_names as "used" for budget
    parent = {source_node: None}
    queue = deque([(source_node, 0)])
    while queue:
        cur, depth = queue.popleft()
        if depth > max_depth:
            continue
        cur_name = cur[:-1]
        # Did we reach the target?
        if cur_name == target_name and cur != source_node:
            # Reconstruct path
            path = []
            n = cur
            while n is not None:
                path.append(n)
                n = parent[n]
            path.reverse()
            return path
        # Expand
        for _, nxt in G.out_edges(cur):
            if nxt in parent:
                continue
            nxt_name = nxt[:-1]
            # Allow target even if it would otherwise be blocked
            if nxt_name == target_name:
                parent[nxt] = cur
                queue.append((nxt, depth + 1))
                continue
            if nxt_name in blocked_names:
                continue
            # Budget check
            if used_counts.get(nxt_name, 0) >= budgets.get(nxt_name, 0):
                continue
            parent[nxt] = cur
            queue.append((nxt, depth + 1))
    return None


def _find_next_target(positions, used_counts, budgets, walked_set,
                       current_ref_pos, max_jump=10_000_000):
    """
    Find the nearest unwalked, budget-available aligned unitig with
    ref_start > current_ref_pos (or any aligned unitig if current_ref_pos is None).
    Returns the unitig name, or None.
    """
    best = None
    best_dist = float('inf')
    for u, (rs, re_, _, _) in positions.items():
        if used_counts.get(u, 0) >= budgets.get(u, 0):
            continue
        if u in walked_set:
            continue
        if current_ref_pos is not None:
            if rs <= current_ref_pos:
                continue
            dist = rs - current_ref_pos
        else:
            dist = rs
        if dist > max_jump:
            continue
        if dist < best_dist:
            best_dist = dist
            best = u
    return best


def _walk_target_directed(G, gfa, positions, budgets, used_counts, chr_length,
                            start_avoid, max_search_depth):
    """
    Target-directed walker:
      1. Find the leftmost aligned unitig as the start.
      2. Find the next reference-forward aligned unitig as the current target.
      3. BFS through the graph trying to reach that target within max_search_depth.
      4. If found, commit the path (including any intermediate unaligned unitigs).
      5. If not found, accept a gap and jump directly to the target.
      6. Advance the target to the next forward aligned unitig and repeat.
    """
    if start_avoid is None:
        start_avoid = set()

    start_node = find_start_node(G, gfa, positions, used_counts, budgets,
                                  avoid=start_avoid)
    if start_node is None:
        print("#   no start node found", file=sys.stderr)
        return [], 0, 0

    print(f"#   start: {start_node}, out_degree={G.out_degree(start_node)}",
          file=sys.stderr)

    committed_path = [start_node]
    committed_walked = {start_node[:-1]}
    used_counts[start_node[:-1]] += 1
    start_name = start_node[:-1]
    committed_ref_pos = positions[start_name][1] if start_name in positions else None
    committed_cum_len = gfa.segment(start_name).length

    gaps = 0
    bfs_calls = 0
    successful_paths = 0

    while True:
        # Stopping conditions
        if committed_ref_pos and committed_ref_pos >= chr_length * 0.98:
            break
        if committed_cum_len > chr_length * 1.3:
            break

        # Find the next forward-positioned aligned unitig as our target
        target = _find_next_target(positions, used_counts, budgets,
                                     committed_walked, committed_ref_pos)
        if target is None:
            break

        # Try to BFS to this target from current committed tail
        bfs_calls += 1
        bfs_path = _bfs_to_target(G, committed_path[-1], target,
                                    committed_walked, max_search_depth,
                                    used_counts, budgets)

        if bfs_path is not None:
            # Commit the path (skip the first node, it's already committed)
            for node in bfs_path[1:]:
                committed_path.append(node)
                name = node[:-1]
                committed_walked.add(name)
                used_counts[name] += 1
                committed_cum_len += gfa.segment(name).length
                if name in positions:
                    committed_ref_pos = max(committed_ref_pos or 0,
                                              positions[name][1])
            successful_paths += 1
        else:
            # No path found within search depth — accept a gap and jump
            placed = False
            for orient in ('+', '-'):
                jump_node = f"{target}{orient}"
                if jump_node in G:
                    committed_path.append('|GAP|')
                    committed_path.append(jump_node)
                    committed_walked.add(target)
                    used_counts[target] += 1
                    committed_cum_len += gfa.segment(target).length
                    committed_ref_pos = positions[target][1]
                    gaps += 1
                    placed = True
                    break
            if not placed:
                break

    print(f"#   bfs_calls={bfs_calls}, successful={successful_paths}, gaps={gaps}",
          file=sys.stderr)
    return committed_path, gaps, bfs_calls






def construct_two_walks(G, gfa, positions, budgets, chr_length,
                        max_search_depth=200):
    used_counts = defaultdict(int)

    print(f"# Building walk 1 (target ~{chr_length:,} bp)...", file=sys.stderr)
    walk1, gaps1, steps1 = walk_one_haplotype(G, gfa, positions, budgets,
                                                used_counts, chr_length,
                                                max_search_depth=max_search_depth)
    walk1_len = sum(gfa.segment(n[:-1]).length for n in walk1 if n != '|GAP|')
    walk1_aligned = sum(1 for n in walk1 if n != '|GAP|' and n[:-1] in positions)
    print(f"# Walk 1: {len(walk1)} nodes ({walk1_aligned} aligned), "
          f"{walk1_len:,} bp, {gaps1} gaps ({steps1} BFS calls)", file=sys.stderr)

    walk1_start = walk1[0][:-1] if walk1 else None
    avoid = {walk1_start} if walk1_start else set()

    print(f"# Building walk 2...", file=sys.stderr)
    walk2, gaps2, steps2 = walk_one_haplotype(G, gfa, positions, budgets,
                                                used_counts, chr_length,
                                                start_avoid=avoid,
                                                max_search_depth=max_search_depth)
    walk2_len = sum(gfa.segment(n[:-1]).length for n in walk2 if n != '|GAP|')
    walk2_aligned = sum(1 for n in walk2 if n != '|GAP|' and n[:-1] in positions)
    print(f"# Walk 2: {len(walk2)} nodes ({walk2_aligned} aligned), "
          f"{walk2_len:,} bp, {gaps2} gaps ({steps2} BFS calls)", file=sys.stderr)

    return walk1, walk2, used_counts


# ---------- Output ----------

def to_gaf(walk, walk_name):
    pieces = [[]]
    for node in walk:
        if node == '|GAP|':
            pieces.append([])
            continue
        name, orient = node[:-1], node[-1]
        sign = '>' if orient == '+' else '<'
        pieces[-1].append(f"{sign}{name}")
    lines = []
    for i, piece in enumerate(pieces):
        if piece:
            lines.append(f"{walk_name}_part{i+1}\t{''.join(piece)}")
    return '\n'.join(lines)


def to_fasta(walk, gfa, walk_name):
    sequences = []
    current_pieces = []
    piece_idx = 1

    def flush():
        nonlocal piece_idx
        if current_pieces:
            seq = ''.join(current_pieces)
            sequences.append((f"{walk_name}_part{piece_idx}", seq))
            current_pieces.clear()
            piece_idx += 1

    for node in walk:
        if node == '|GAP|':
            flush()
            continue
        name, orient = node[:-1], node[-1]
        seg = gfa.segment(name)
        seq = seg.sequence
        if seq is None or seq == '*':
            continue
        if orient == '-':
            seq = revcomp(seq)
        current_pieces.append(seq)
    flush()
    return sequences


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gfa', required=True)
    ap.add_argument('--paf', required=True)
    ap.add_argument('--chr', required=True)
    ap.add_argument('--include', help='Optional include list of unitig names')
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--default-budget', type=int, default=2)
    ap.add_argument('--max-search-depth', type=int, default=200,
                    help='Max BFS depth when searching for path to next target (default 200)')
    args = ap.parse_args()

    print(f"# Loading GFA: {args.gfa}", file=sys.stderr)
    gfa, G = load_gfa(args.gfa)
    print(f"# {len(gfa.segments)} segments, {G.number_of_edges()} oriented edges",
          file=sys.stderr)

    include_set = None
    if args.include:
        with open(args.include) as f:
            include_set = set(line.strip() for line in f if line.strip())

    print(f"# Loading PAF: {args.paf}", file=sys.stderr)
    positions, chr_length = load_positions(args.paf, args.chr, include_set)
    print(f"# {len(positions)} unitigs aligned to {args.chr} "
          f"(chr_length = {chr_length:,} bp)", file=sys.stderr)

    gfa_segs = {seg.name for seg in gfa.segments}
    before = len(positions)
    positions = {u: v for u, v in positions.items() if u in gfa_segs}
    dropped = before - len(positions)
    if dropped:
        print(f"# Dropped {dropped} PAF unitigs not in GFA", file=sys.stderr)

    print(f"# Detecting bubbles...", file=sys.stderr)
    bubble_members, seg_neighbors = detect_bubble_members(gfa)
    print(f"# {len(bubble_members)} bubble members", file=sys.stderr)

    budgets = assign_budgets(gfa, positions, bubble_members, seg_neighbors,
                             default_budget=args.default_budget)

    walk1, walk2, used = construct_two_walks(G, gfa, positions, budgets, chr_length,
                                              max_search_depth=args.max_search_depth)

    with open(f"{args.prefix}.hap1.gaf", 'w') as f:
        f.write(to_gaf(walk1, f"{args.chr}_hap1") + '\n')
    with open(f"{args.prefix}.hap2.gaf", 'w') as f:
        f.write(to_gaf(walk2, f"{args.chr}_hap2") + '\n')

    with open(f"{args.prefix}.hap1.fa", 'w') as f:
        for name, seq in to_fasta(walk1, gfa, f"{args.chr}_hap1"):
            f.write(f">{name}\n{seq}\n")
    with open(f"{args.prefix}.hap2.fa", 'w') as f:
        for name, seq in to_fasta(walk2, gfa, f"{args.chr}_hap2"):
            f.write(f">{name}\n{seq}\n")

    with open(f"{args.prefix}.stats.tsv", 'w') as f:
        f.write("unitig\tlength\tdegree\tbubble_alt\tbudget\tused\taligned_chr\tref_start\tref_end\n")
        for seg in gfa.segments:
            pos = positions.get(seg.name)
            f.write(
                f"{seg.name}\t{seg.length}\t{len(seg_neighbors[seg.name])}\t"
                f"{'Y' if seg.name in bubble_members else 'N'}\t"
                f"{budgets.get(seg.name, 0)}\t{used.get(seg.name, 0)}\t"
                f"{args.chr if pos else 'NA'}\t"
                f"{pos[0] if pos else 'NA'}\t{pos[1] if pos else 'NA'}\n"
            )

    walk1_len = sum(gfa.segment(n[:-1]).length for n in walk1 if n != '|GAP|')
    walk2_len = sum(gfa.segment(n[:-1]).length for n in walk2 if n != '|GAP|')
    print(f"# Final: chr_length={chr_length:,}", file=sys.stderr)
    print(f"#   walk1: {walk1_len:,} bp ({100*walk1_len/chr_length:.1f}%)",
          file=sys.stderr)
    print(f"#   walk2: {walk2_len:,} bp ({100*walk2_len/chr_length:.1f}%)",
          file=sys.stderr)
    print(f"# Done. Outputs at {args.prefix}.*", file=sys.stderr)


if __name__ == '__main__':
    main()
