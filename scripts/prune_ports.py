#!/usr/bin/env python3
"""
prune_ports.py
--------------
Enforce the topology rule that each node port has at most one incident
edge, EXCEPT where the multiple edges form a bubble that reconverges
somewhere downstream.

Bidirected ports:
    Every node has two ports: L (left side) and R (right side).
    An L line "A oA B oB" attaches to:
        - port_R of A if oA == '+', else port_L of A
        - port_L of B if ob == '+', else port_R of B

For each port with >1 edges, the script computes the forward-reachable
set of (node, entered_port) states from each edge's target. If all pairs
of edges at this port have non-empty intersections, the multi-edge
configuration represents a legitimate bubble (alternative paths
reconverge in the graph), and ALL edges are kept.

If any pair of edges has empty intersection (alternative paths never
reconverge), the script declares this a structural break and cuts edges
to reduce the port to one. The edge KEPT at each such port is the one
whose downstream chain has the greatest greedy length, capped at
MAX_WALK_NODES.

This is purely topological -- no length or coverage thresholds determine
whether a multi-edge port is a bubble or a break.

Usage:
    python prune_ports.py <in.gfa> <out.gfa> [--max-walk-nodes N] [--verbose]
"""

import argparse
import sys
from collections import defaultdict, deque


MAX_WALK_NODES_DEFAULT = 100


def parse_gfa(path):
    segments = {}
    links = []
    raw_lines = []
    with open(path) as fh:
        for line in fh:
            raw_lines.append(line)
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                segments[parts[1]] = parts[2]
            elif line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                a, oa, b, ob = parts[1], parts[2], parts[3], parts[4]
                rest = parts[5:] if len(parts) > 5 else ["0M"]
                links.append((a, oa, b, ob, rest))
    return segments, links, raw_lines


def flip(o):
    return "-" if o == "+" else "+"


def edge_endpoints(a, oa, b, ob):
    a_port = "R" if oa == "+" else "L"
    b_port = "L" if ob == "+" else "R"
    return (a, a_port), (b, b_port)


def canonical_edge_key(a, oa, b, ob):
    forward = (a, oa, b, ob)
    reverse = (b, flip(ob), a, flip(oa))
    return min(forward, reverse)


def build_port_index(links):
    port_edges = defaultdict(list)
    seen_edges = set()
    for a, oa, b, ob, _rest in links:
        key = canonical_edge_key(a, oa, b, ob)
        if key in seen_edges:
            continue
        seen_edges.add(key)
        e1, e2 = edge_endpoints(a, oa, b, ob)
        port_edges[e1].append(key)
        port_edges[e2].append(key)
    return port_edges


def build_adjacency(links):
    """
    adj[(node, port)] = list of (neighbor_node, neighbor_entered_port, edge_key)
    where neighbor_entered_port is the port we ARRIVE at on the neighbor.
    """
    adj = defaultdict(list)
    seen_edges = set()
    for a, oa, b, ob, _rest in links:
        key = canonical_edge_key(a, oa, b, ob)
        if key in seen_edges:
            continue
        seen_edges.add(key)
        (na, pa), (nb, pb) = edge_endpoints(a, oa, b, ob)
        adj[(na, pa)].append((nb, pb, key))
        adj[(nb, pb)].append((na, pa, key))
    return adj


def forward_reachable(start_node, entered_port, adj, cache):
    """
    BFS forward from a state where we have just entered start_node at
    entered_port. We then exit via the OPPOSITE port and continue.

    Returns the set of (node, entered_port) states reachable, including
    the starting state itself.

    Memoized in `cache` keyed by (start_node, entered_port).
    """
    key = (start_node, entered_port)
    if key in cache:
        return cache[key]

    reachable = set()
    queue = deque([(start_node, entered_port)])
    reachable.add(key)

    while queue:
        node, enter_port = queue.popleft()
        exit_port = "L" if enter_port == "R" else "R"
        for nb_node, nb_enter, _ek in adj.get((node, exit_port), []):
            state = (nb_node, nb_enter)
            if state not in reachable:
                reachable.add(state)
                queue.append(state)

    cache[key] = reachable
    return reachable


def chain_length_from(start_node, entry_port, segments, adj, max_nodes):
    """
    Greedy walk used only for tiebreaking when we have to cut.
    Same semantics as forward_reachable but greedy and length-aware.
    """
    visited = set()
    total = 0
    node = start_node
    enter_port = entry_port

    for _ in range(max_nodes):
        if node in visited:
            break
        visited.add(node)
        total += len(segments.get(node, ""))

        exit_port = "L" if enter_port == "R" else "R"
        options = adj.get((node, exit_port), [])
        if not options:
            break

        def opt_key(opt):
            nb_node, nb_enter, _ = opt
            return (-len(segments.get(nb_node, "")), nb_node, nb_enter)

        best = min(options, key=opt_key)
        node, enter_port, _ = best

    return total


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("gfa_in")
    ap.add_argument("gfa_out")
    ap.add_argument("--max-walk-nodes", type=int, default=MAX_WALK_NODES_DEFAULT,
                    help=f"Cap greedy chain walk for tiebreaking at N nodes "
                         f"(default {MAX_WALK_NODES_DEFAULT})")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    segments, links, raw_lines = parse_gfa(args.gfa_in)
    if args.verbose:
        print(f"[load] {len(segments)} segments, {len(links)} links",
              file=sys.stderr)

    links_by_key = {}
    for link in links:
        a, oa, b, ob, _ = link
        key = canonical_edge_key(a, oa, b, ob)
        links_by_key[key] = link

    port_edges = build_port_index(links)
    adj = build_adjacency(links)

    over_ports = [(p, edges) for p, edges in port_edges.items() if len(edges) > 1]
    if args.verbose:
        print(f"[scan] {len(over_ports)} ports have >1 edges", file=sys.stderr)

    reach_cache = {}
    edges_to_cut = set()
    n_kept_as_bubble = 0
    n_cut_as_break = 0

    for (port_node, port_side), edge_keys in over_ports:
        # For each edge at this port, find the (other_node, other_port)
        # we'd be entering. That's the BFS starting state.
        edge_targets = []
        for key in edge_keys:
            a, oa, b, ob, _ = links_by_key[key]
            (na, pa), (nb, pb) = edge_endpoints(a, oa, b, ob)
            if (na, pa) == (port_node, port_side):
                target = (nb, pb)
            else:
                target = (na, pa)
            edge_targets.append((key, target))

        # Compute reachable sets per edge.
        reach_sets = {}
        for key, (tn, tp) in edge_targets:
            reach_sets[key] = forward_reachable(tn, tp, adj, reach_cache)

        # Build the "reconverges with" graph among edges at this port:
        # nodes are edges; an undirected edge connects two edges whose
        # reachable sets intersect. Then find connected components.
        from collections import defaultdict as _dd
        groups = {key: key for key, _ in edge_targets}  # union-find

        def find(x):
            while groups[x] != x:
                groups[x] = groups[groups[x]]
                x = groups[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                groups[rx] = ry

        keys_list = [k for k, _ in edge_targets]
        for i in range(len(keys_list)):
            for j in range(i + 1, len(keys_list)):
                if reach_sets[keys_list[i]] & reach_sets[keys_list[j]]:
                    union(keys_list[i], keys_list[j])

        # Group by root.
        clusters = _dd(list)
        for k in keys_list:
            clusters[find(k)].append(k)

        # If there's exactly one cluster and it contains all edges, this
        # is a clean bubble -- keep all edges.
        if len(clusters) == 1 and len(next(iter(clusters.values()))) == len(keys_list):
            n_kept_as_bubble += 1
            continue

        # Otherwise we have at least one edge that doesn't reconverge
        # with all the others. Keep the LARGEST cluster (the cluster
        # containing the most edges -- representing the largest coherent
        # bubble or chain at this port). Cut all edges not in it.
        # Tiebreak between equally-sized clusters: pick the one whose
        # edges' targets give the largest total greedy chain length.
        def cluster_score(edge_list):
            total = 0
            for key in edge_list:
                tn, tp = dict(edge_targets)[key]
                total += chain_length_from(tn, tp, segments, adj,
                                            args.max_walk_nodes)
            return (len(edge_list), total)

        best_cluster = max(clusters.values(), key=cluster_score)

        n_cut_as_break += 1
        for k in keys_list:
            if k not in best_cluster:
                edges_to_cut.add(k)

    if args.verbose:
        print(f"[decide] {n_kept_as_bubble} multi-edge ports preserved "
              "as bubbles", file=sys.stderr)
        print(f"[decide] {n_cut_as_break} multi-edge ports cut as "
              "structural breaks", file=sys.stderr)
        print(f"[decide] {len(edges_to_cut)} edges marked for cutting",
              file=sys.stderr)

    # Write the GFA.
    n_kept = n_cut = 0
    with open(args.gfa_out, "w") as fout:
        for line in raw_lines:
            if not line.startswith("L\t"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            a, oa, b, ob = parts[1], parts[2], parts[3], parts[4]
            key = canonical_edge_key(a, oa, b, ob)
            if key in edges_to_cut:
                n_cut += 1
            else:
                fout.write(line)
                n_kept += 1

    print(f"[done] kept {n_kept} links, cut {n_cut} -> {args.gfa_out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
