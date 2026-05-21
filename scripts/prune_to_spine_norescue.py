#!/usr/bin/env python3
"""
prune_to_spine.py
-----------------
For each connected component of a GFA, identify the two "true tips" (the
pair of degree-1 nodes whose base-pair-weighted shortest path between
them is maximal -- i.e., the two endpoints of the component's diameter
measured in bases) and remove every node not on at least one simple
path between them.

The spine endpoints are chosen by computing single-source weighted
shortest paths from every degree-1 node (Dijkstra) and selecting the
pair maximizing the distance. This finds the real chromosome tips even
when the graph has internal tangles with their own dead-end tips: an
internal dead-end is graph-close to all other tips because the tangle
is locally connected, whereas the real chromosome ends are graph-far
because they're separated by the chromosome's length.

Bubbles are preserved -- nodes on either side of a bubble between the
spine tips lie on simple paths between them and survive the BCC step.

Per-component algorithm:
    1. Build undirected graph of the component.
    2. Find all degree-1 nodes; these are tip candidates.
    3. For each tip, run Dijkstra to all other nodes (weight = node length
       added at arrival). Among (tip, tip) pairs, pick the one with the
       greatest distance: this is the spine.
    4. Use biconnected components (with a virtual edge between t1* and
       t2*) to identify all nodes lying on SOME simple path between
       them -- this is the KEEP set.
    5. Drop every other node.

Components with < 2 degree-1 nodes (pure cycles, isolated nodes) are
passed through unchanged with a warning.

Dependencies: networkx.

Usage:
    python prune_to_spine.py <in.gfa> <out.gfa> \\
        [--tip-tolerance BP] [--min-component-size N] [--verbose]
"""

import argparse
import sys
from collections import defaultdict

try:
    import networkx as nx
except ImportError:
    sys.exit("ERROR: this script requires networkx. Install with:\n"
             "    pip install networkx")


def parse_gfa(path):
    """
    Return (segments, links, raw_lines).
        segments: dict node_id -> sequence
        links: list of (a, oa, b, ob, cigar, tags) tuples
        raw_lines: original line list (preserves order for write-back)
    """
    segments = {}
    links = []
    raw_lines = []
    with open(path) as fh:
        for line in fh:
            raw_lines.append(line)
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                nid, seq = parts[1], parts[2]
                segments[nid] = seq
            elif line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                a, oa, b, ob = parts[1], parts[2], parts[3], parts[4]
                cigar = parts[5] if len(parts) > 5 else "0M"
                tags = parts[6:] if len(parts) > 6 else []
                links.append((a, oa, b, ob, cigar, tags))
    return segments, links, raw_lines


def build_undirected_graph(segments, links):
    """Build an undirected networkx graph; edges weighted by neighbor length."""
    G = nx.Graph()
    for nid, seq in segments.items():
        G.add_node(nid, length=len(seq))
    for a, _, b, _, _, _ in links:
        if a in segments and b in segments:
            # Multi-edges (rare in GFA) get collapsed to a single edge here;
            # that's fine for spine and BCC analysis.
            G.add_edge(a, b)
    return G


def degree_one_nodes(G):
    """Undirected degree-1 nodes. Kept for backwards compat / fallback use."""
    return [n for n in G.nodes() if G.degree(n) == 1]


def find_bidirected_tips(segments, links):
    """
    A node is a "tip" in the bidirected sense if at least one of its two
    ports (left or right) has zero incident edges. Equivalently, all edges
    incident to the node attach to the same port -- the chromosome ends
    there, even if there are multiple edges leaving the same side.

    GFA port mapping for link "A oA B oB":
        - port_R of A if oA == '+', else port_L of A
        - port_L of B if ob == '+', else port_R of B
    """
    port_count = defaultdict(lambda: {"L": 0, "R": 0})
    for a, oa, b, ob, *_ in links:
        if a not in segments or b not in segments:
            continue
        a_port = "R" if oa == "+" else "L"
        b_port = "L" if ob == "+" else "R"
        port_count[a][a_port] += 1
        port_count[b][b_port] += 1

    tips = []
    for nid in segments:
        ports = port_count[nid]
        # Tip if at least one port has zero edges (and the node has any
        # edges at all -- pure isolated nodes are handled separately).
        if (ports["L"] == 0) != (ports["R"] == 0):
            tips.append(nid)
        elif ports["L"] == 0 and ports["R"] == 0:
            # Isolated node: both ports empty. Skip -- not a useful spine
            # endpoint.
            pass
    return tips


def find_spine(G, tips, total_time_budget=None, verbose=False):
    """
    Among all degree-1 nodes, find the pair (t1, t2) maximizing the
    base-pair-weighted shortest path between them in the graph. That pair
    is the "graph diameter" measured in bases.

    Returns (t1, t2, length, second_best_length).

    This is exact, fast, and finds the chromosome spine even when the
    graph has tangled internal structure. The intuition: a tip that
    represents a real chromosome end is graph-far from the other real
    end (separated by the entire chromosome length); a spurious tip
    inside a tangle is graph-close to all other tips because tangles
    are locally connected.

    total_time_budget is ignored (the algorithm is so fast it's not
    needed); kept in the signature for backwards compatibility.
    """
    if len(tips) < 2:
        sys.exit("ERROR: graph has fewer than 2 degree-1 nodes; "
                 "cannot identify a tip pair.")

    tips_set = set(tips)
    pair_lengths = {}  # frozenset({a, b}) -> distance

    # Single-source weighted shortest paths from each tip, using node
    # length as the weight added when ARRIVING at the node. We use a
    # weight function that returns the length of node v (the target).
    def weight_fn(u, v, _data):
        return G.nodes[v]["length"]

    for t in tips:
        if verbose:
            print(f"[spine] Dijkstra from tip {t}", file=sys.stderr)
        distances = nx.single_source_dijkstra_path_length(
            G, t, weight=weight_fn
        )
        # Add source's own length to get full path length including both endpoints.
        src_len = G.nodes[t]["length"]
        for v, d in distances.items():
            if v in tips_set and v != t:
                full_len = src_len + d
                pair = frozenset({t, v})
                if full_len > pair_lengths.get(pair, 0):
                    pair_lengths[pair] = full_len

    if not pair_lengths:
        sys.exit("ERROR: no tip-to-tip path exists in the component.")

    ranked = sorted(pair_lengths.items(), key=lambda kv: -kv[1])
    best_pair, best_len = ranked[0]
    second_best_len = ranked[1][1] if len(ranked) > 1 else 0
    t1, t2 = sorted(best_pair)

    return t1, t2, best_len, second_best_len


def nodes_on_paths_between(G, t1, t2):
    """
    Return the set of nodes lying on at least one simple path between
    t1 and t2 in G.

    Method: add a virtual edge (t1, t2), compute biconnected components,
    return all nodes in the BCC containing the virtual edge.

    This works because two nodes lie in the same biconnected component
    iff there are two internally-vertex-disjoint paths between them
    (Menger's theorem). With the virtual edge present, the BCC containing
    it is exactly the set of nodes that can reach both t1 and t2 along
    internally-disjoint paths -- equivalently, that lie on some simple
    path between t1 and t2.
    """
    H = G.copy()
    H.add_edge(t1, t2, virtual=True)
    for bcc in nx.biconnected_components(H):
        if t1 in bcc and t2 in bcc:
            return set(bcc)
    # Fallback should never happen if t1, t2 are connected.
    return {t1, t2}


def write_filtered_gfa(raw_lines, keep_nodes, out_path, verbose=False):
    """Write a new GFA dropping S lines and L lines that reference dropped nodes."""
    n_seg_kept = 0
    n_seg_dropped = 0
    n_link_kept = 0
    n_link_dropped = 0

    with open(out_path, "w") as fout:
        for line in raw_lines:
            if line.startswith("S\t"):
                nid = line.split("\t", 2)[1]
                if nid in keep_nodes:
                    fout.write(line)
                    n_seg_kept += 1
                else:
                    n_seg_dropped += 1
            elif line.startswith("L\t"):
                parts = line.split("\t")
                a, b = parts[1], parts[3]
                if a in keep_nodes and b in keep_nodes:
                    fout.write(line)
                    n_link_kept += 1
                else:
                    n_link_dropped += 1
            else:
                fout.write(line)

    if verbose:
        print(f"[write] kept {n_seg_kept} segments, dropped {n_seg_dropped}",
              file=sys.stderr)
        print(f"[write] kept {n_link_kept} links, dropped {n_link_dropped}",
              file=sys.stderr)


def process_component(G, comp_id, tips_in_component, args):
    """
    Process one component and return its keep-set (a set of node IDs).

    tips_in_component is the list of bidirected tips that belong to this
    component (a subset of the global tip list).

    Returns (keep_set, summary_dict) where summary_dict has keys:
        size, n_tips, t1, t2, spine_len, second_best, ambiguous, status
    status is one of: "ok", "too_few_tips", "too_small"
    """
    n_nodes = G.number_of_nodes()
    summary = {"comp_id": comp_id, "size": n_nodes,
               "t1": None, "t2": None, "spine_len": 0,
               "second_best": 0, "ambiguous": False, "status": "ok",
               "n_tips": len(tips_in_component)}

    if n_nodes < args.min_component_size:
        summary["status"] = "too_small"
        return set(G.nodes()), summary

    if len(tips_in_component) < 2:
        # No spine possible (pure cycle, isolated, or only one tip).
        summary["status"] = "too_few_tips"
        return set(G.nodes()), summary

    t1, t2, spine_len, second_best = find_spine(
        G, tips_in_component, verbose=args.verbose,
    )

    summary["t1"] = t1
    summary["t2"] = t2
    summary["spine_len"] = spine_len
    summary["second_best"] = second_best
    summary["ambiguous"] = (
        second_best > 0
        and spine_len - second_best < args.tip_tolerance
    )

    keep = nodes_on_paths_between(G, t1, t2)
    return keep, summary


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("gfa_in")
    ap.add_argument("gfa_out")
    ap.add_argument("--tip-tolerance", type=int, default=10000,
                    help="If multiple tip pairs are within TIP_TOLERANCE bp "
                         "of the longest, warn (in bp). Default: 10000.")
    ap.add_argument("--max-search-time", type=float, default=60.0,
                    help="Ignored (kept for backwards compatibility). "
                         "The spine search is now exact and runs in a few "
                         "seconds even on large graphs.")
    ap.add_argument("--min-component-size", type=int, default=1,
                    help="Skip pruning components smaller than this many "
                         "nodes; pass them through unchanged. Default: 1.")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    segments, links, raw_lines = parse_gfa(args.gfa_in)
    G_all = build_undirected_graph(segments, links)
    print(f"[load] {len(segments)} segments, {G_all.number_of_edges()} edges",
          file=sys.stderr)

    # Identify bidirected tips: nodes with an empty port (left or right).
    # This is more general than undirected degree-1, which misses nodes
    # whose edges all attach to the same side (those are also tips since
    # the chromosome ends there).
    all_tips = set(find_bidirected_tips(segments, links))
    print(f"[load] {len(all_tips)} bidirected tips", file=sys.stderr)

    components = sorted(nx.connected_components(G_all), key=len, reverse=True)
    print(f"[load] {len(components)} connected component(s)", file=sys.stderr)

    # Process each component independently. Union all keep-sets.
    union_keep = set()
    summaries = []
    for i, comp_nodes in enumerate(components):
        sub = G_all.subgraph(comp_nodes).copy()
        comp_tips = [t for t in all_tips if t in comp_nodes]
        keep, summary = process_component(sub, i, comp_tips, args)
        union_keep |= keep
        summaries.append(summary)

        if summary["status"] == "ok":
            note = " (ambiguous)" if summary["ambiguous"] else ""
            print(f"[comp {i}] size={summary['size']} tips={summary['n_tips']} "
                  f"spine={summary['t1']}<->{summary['t2']} "
                  f"len={summary['spine_len']}bp "
                  f"kept={len(keep)}{note}", file=sys.stderr)
        elif summary["status"] == "too_few_tips":
            print(f"[comp {i}] size={summary['size']} "
                  f"tips={summary['n_tips']} -- "
                  "fewer than 2 tips, passing through unmodified",
                  file=sys.stderr)
        elif summary["status"] == "too_small":
            print(f"[comp {i}] size={summary['size']} "
                  "below --min-component-size, passing through unmodified",
                  file=sys.stderr)

    write_filtered_gfa(raw_lines, union_keep, args.gfa_out,
                       verbose=args.verbose)

    # Write the tips file: <gfa_out>.tips
    # Format: TSV with columns
    #   comp_id  status  size  n_tips  tip1  tip2  spine_len_bp  ambiguous  second_best_bp
    # status is "ok" / "too_few_tips" / "too_small".
    tips_path = args.gfa_out + ".tips"
    with open(tips_path, "w") as tfh:
        tfh.write("\t".join([
            "comp_id", "status", "size", "n_tips",
            "tip1", "tip2", "spine_len_bp",
            "ambiguous", "second_best_bp",
        ]) + "\n")
        for s in summaries:
            tfh.write("\t".join([
                str(s["comp_id"]),
                s["status"],
                str(s["size"]),
                str(s["n_tips"]),
                str(s["t1"]) if s["t1"] is not None else ".",
                str(s["t2"]) if s["t2"] is not None else ".",
                str(s["spine_len"]),
                "yes" if s["ambiguous"] else "no",
                str(s["second_best"]),
            ]) + "\n")
    print(f"[done] wrote per-component tips to {tips_path}", file=sys.stderr)

    # Summary stats.
    total_in = len(segments)
    total_kept = len(union_keep)
    n_ok = sum(1 for s in summaries if s["status"] == "ok")
    n_ambig = sum(1 for s in summaries if s["ambiguous"])
    print(f"[done] {n_ok}/{len(components)} components had a spine identified "
          f"({n_ambig} ambiguous)", file=sys.stderr)
    print(f"[done] kept {total_kept} / {total_in} segments "
          f"({100.0 * total_kept / max(total_in, 1):.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
