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

Rescue passes:
    After the initial round, the nodes that were dropped may themselves
    contain large tip-to-tip components (e.g. a second chromosome arm
    that was attached to the main spine through a tangle). With
    --rescue-min-length set, the script will re-run the spine-finding
    algorithm on the leftover subgraph and keep any newly-found spines
    whose length meets the threshold. Rescue runs iteratively until no
    component meets the threshold or --max-rescue-rounds is reached.
    Links between the original spine and rescued components are dropped
    naturally (any L-line referencing a still-dropped node is filtered),
    so the rescued components come out as separate connected pieces in
    the output GFA -- ready for independent assembly.

Dependencies: networkx.

Usage:
    python prune_to_spine.py <in.gfa> <out.gfa> \\
        [--tip-tolerance BP] [--min-component-size N] \\
        [--rescue-min-length BP] [--max-rescue-rounds N] [--verbose]
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


def build_undirected_graph(segments, links, node_subset=None,
                           link_subset=None):
    """
    Build an undirected networkx graph.

    If node_subset is provided, only those nodes are added.
    If link_subset is provided, only those links are considered; otherwise
    all links among included nodes are added.
    """
    G = nx.Graph()
    nodes_to_use = node_subset if node_subset is not None else segments.keys()
    for nid in nodes_to_use:
        if nid in segments:
            G.add_node(nid, length=len(segments[nid]))

    links_to_use = link_subset if link_subset is not None else links
    node_set = set(G.nodes())
    for a, _, b, _, _, _ in links_to_use:
        if a in node_set and b in node_set:
            G.add_edge(a, b)
    return G


def degree_one_nodes(G):
    """Undirected degree-1 nodes. Kept for backwards compat / fallback use."""
    return [n for n in G.nodes() if G.degree(n) == 1]


def find_bidirected_tips(segments, links, node_subset=None):
    """
    A node is a "tip" in the bidirected sense if at least one of its two
    ports (left or right) has zero incident edges. Equivalently, all edges
    incident to the node attach to the same port -- the chromosome ends
    there, even if there are multiple edges leaving the same side.

    GFA port mapping for link "A oA B oB":
        - port_R of A if oA == '+', else port_L of A
        - port_L of B if ob == '+', else port_R of B

    If node_subset is given, only edges where BOTH endpoints are in the
    subset are counted -- so a node's tip status is evaluated as it
    appears in the induced subgraph. This is what rescue passes need:
    nodes that were internal in the full graph but become tips once
    their neighbors in the kept set have been removed.
    """
    port_count = defaultdict(lambda: {"L": 0, "R": 0})
    in_subset = (lambda nid: nid in node_subset) if node_subset is not None \
        else (lambda nid: nid in segments)

    for a, oa, b, ob, *_ in links:
        if not in_subset(a) or not in_subset(b):
            continue
        a_port = "R" if oa == "+" else "L"
        b_port = "L" if ob == "+" else "R"
        port_count[a][a_port] += 1
        port_count[b][b_port] += 1

    tips = []
    candidate_nodes = node_subset if node_subset is not None else segments.keys()
    for nid in candidate_nodes:
        if nid not in segments:
            continue
        ports = port_count[nid]
        # Tip if exactly one port is empty (one-sided dead end).
        # If both are empty, the node is isolated in this view -- skip;
        # it's not a useful spine endpoint.
        if (ports["L"] == 0) != (ports["R"] == 0):
            tips.append(nid)
    return tips


def find_spine(G, tips, verbose=False):
    """
    Among all tips, find the pair (t1, t2) maximizing the base-pair-
    weighted shortest path between them in the graph.

    Returns (t1, t2, length, second_best_length).
    """
    if len(tips) < 2:
        return None, None, 0, 0

    tips_set = set(tips)
    pair_lengths = {}  # frozenset({a, b}) -> distance

    def weight_fn(u, v, _data):
        return G.nodes[v]["length"]

    for t in tips:
        if verbose:
            print(f"[spine] Dijkstra from tip {t}", file=sys.stderr)
        distances = nx.single_source_dijkstra_path_length(
            G, t, weight=weight_fn
        )
        src_len = G.nodes[t]["length"]
        for v, d in distances.items():
            if v in tips_set and v != t:
                full_len = src_len + d
                pair = frozenset({t, v})
                if full_len > pair_lengths.get(pair, 0):
                    pair_lengths[pair] = full_len

    if not pair_lengths:
        return None, None, 0, 0

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
    """
    H = G.copy()
    H.add_edge(t1, t2, virtual=True)
    for bcc in nx.biconnected_components(H):
        if t1 in bcc and t2 in bcc:
            return set(bcc)
    return {t1, t2}


def write_filtered_gfa(raw_lines, keep_nodes, node_group, out_path,
                       verbose=False):
    """
    Write a new GFA dropping S lines and L lines that reference dropped nodes.

    node_group: dict node_id -> group label (typically a comp_id string like
    "0.0" or "1.0"). Links whose two endpoints belong to DIFFERENT groups
    are dropped, so rescued components come out as separate connected
    pieces in the output.

    If a kept node doesn't appear in node_group (defensive fallback), it
    gets group None and any links to it are kept (same-group test still
    succeeds against another None).
    """
    n_seg_kept = 0
    n_seg_dropped = 0
    n_link_kept = 0
    n_link_dropped = 0
    n_link_cross_group = 0

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
                    ga = node_group.get(a)
                    gb = node_group.get(b)
                    if ga == gb:
                        fout.write(line)
                        n_link_kept += 1
                    else:
                        n_link_dropped += 1
                        n_link_cross_group += 1
                else:
                    n_link_dropped += 1
            else:
                fout.write(line)

    if verbose:
        print(f"[write] kept {n_seg_kept} segments, dropped {n_seg_dropped}",
              file=sys.stderr)
        print(f"[write] kept {n_link_kept} links, dropped {n_link_dropped} "
              f"({n_link_cross_group} severed between spine groups)",
              file=sys.stderr)


def process_component(G, comp_id, tips_in_component, args, round_id=0):
    """
    Process one component and return its keep-set (a set of node IDs).

    Returns (keep_set, summary_dict).
    """
    n_nodes = G.number_of_nodes()
    summary = {"comp_id": comp_id, "round": round_id, "size": n_nodes,
               "t1": None, "t2": None, "spine_len": 0,
               "second_best": 0, "ambiguous": False, "status": "ok",
               "n_tips": len(tips_in_component)}

    if n_nodes < args.min_component_size:
        summary["status"] = "too_small"
        return set(G.nodes()) if round_id == 0 else set(), summary

    if len(tips_in_component) < 2:
        summary["status"] = "too_few_tips"
        # Round 0: pass through unchanged (preserve original behavior).
        # Rescue rounds: drop -- a component with <2 tips has no spine,
        # so there's nothing to rescue.
        return (set(G.nodes()) if round_id == 0 else set()), summary

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

    # Rescue threshold: only keep rescued spines that are long enough.
    if round_id > 0 and spine_len < args.rescue_min_length:
        summary["status"] = "below_rescue_threshold"
        return set(), summary

    keep = nodes_on_paths_between(G, t1, t2)
    return keep, summary


def run_pass(segments, links, candidate_nodes, args, round_id):
    """
    Run one full pass of the spine-finding algorithm over every
    connected component of the subgraph induced by candidate_nodes.

    Returns (keep_set, node_group, list_of_summaries) where node_group
    maps each kept node id -> the comp_id of the component it was kept
    by (used downstream to sever links between different spine groups).
    """
    G_all = build_undirected_graph(segments, links, node_subset=candidate_nodes)
    if args.verbose or round_id > 0:
        print(f"[round {round_id}] graph: {G_all.number_of_nodes()} segments, "
              f"{G_all.number_of_edges()} edges", file=sys.stderr)

    all_tips = set(find_bidirected_tips(segments, links, node_subset=candidate_nodes))
    if args.verbose or round_id > 0:
        print(f"[round {round_id}] {len(all_tips)} bidirected tips",
              file=sys.stderr)

    components = sorted(nx.connected_components(G_all), key=len, reverse=True)
    if args.verbose or round_id > 0:
        print(f"[round {round_id}] {len(components)} connected component(s)",
              file=sys.stderr)

    union_keep = set()
    node_group = {}
    summaries = []
    for i, comp_nodes in enumerate(components):
        sub = G_all.subgraph(comp_nodes).copy()
        comp_tips = [t for t in all_tips if t in comp_nodes]
        comp_label = f"{round_id}.{i}"
        keep, summary = process_component(sub, comp_label, comp_tips, args,
                                          round_id=round_id)
        union_keep |= keep
        for nid in keep:
            node_group[nid] = comp_label
        summaries.append(summary)

        cid = summary["comp_id"]
        if summary["status"] == "ok":
            note = " (ambiguous)" if summary["ambiguous"] else ""
            prefix = "[rescue]" if round_id > 0 else "[comp]"
            print(f"{prefix} {cid} size={summary['size']} "
                  f"tips={summary['n_tips']} "
                  f"spine={summary['t1']}<->{summary['t2']} "
                  f"len={summary['spine_len']}bp "
                  f"kept={len(keep)}{note}", file=sys.stderr)
        elif summary["status"] == "too_few_tips":
            if round_id == 0:
                print(f"[comp] {cid} size={summary['size']} "
                      f"tips={summary['n_tips']} -- "
                      "fewer than 2 tips, passing through unmodified",
                      file=sys.stderr)
            elif args.verbose:
                print(f"[rescue] {cid} size={summary['size']} "
                      f"tips={summary['n_tips']} -- "
                      "fewer than 2 tips, skipping",
                      file=sys.stderr)
        elif summary["status"] == "too_small":
            if round_id == 0:
                print(f"[comp] {cid} size={summary['size']} "
                      "below --min-component-size, passing through unmodified",
                      file=sys.stderr)
        elif summary["status"] == "below_rescue_threshold":
            if args.verbose:
                print(f"[rescue] {cid} size={summary['size']} "
                      f"spine={summary['spine_len']}bp < "
                      f"{args.rescue_min_length}bp threshold, skipping",
                      file=sys.stderr)

    return union_keep, node_group, summaries


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
                    help="Ignored (kept for backwards compatibility).")
    ap.add_argument("--min-component-size", type=int, default=1,
                    help="Skip pruning components smaller than this many "
                         "nodes; pass them through unchanged. Default: 1.")
    ap.add_argument("--rescue-min-length", type=int, default=100000,
                    help="After the initial pass, re-run the algorithm on "
                         "dropped nodes and rescue any tip-to-tip component "
                         "with spine length >= this many bp. Set to 0 to "
                         "disable rescue. Default: 100000.")
    ap.add_argument("--max-rescue-rounds", type=int, default=10,
                    help="Maximum number of iterative rescue passes. "
                         "Default: 10. Rescue stops early if a round "
                         "rescues nothing.")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    segments, links, raw_lines = parse_gfa(args.gfa_in)
    print(f"[load] {len(segments)} segments, {len(links)} links",
          file=sys.stderr)

    # ---- Round 0: original pass over all segments. ----
    candidates = set(segments.keys())
    union_keep, node_group, all_summaries = run_pass(
        segments, links, candidates, args, round_id=0)
    round0_kept = len(union_keep)
    print(f"[round 0] kept {round0_kept} / {len(segments)} segments",
          file=sys.stderr)

    # ---- Rescue rounds. ----
    if args.rescue_min_length > 0:
        for round_id in range(1, args.max_rescue_rounds + 1):
            leftover = set(segments.keys()) - union_keep
            if not leftover:
                if args.verbose:
                    print(f"[round {round_id}] no leftover nodes, stopping",
                          file=sys.stderr)
                break

            print(f"[round {round_id}] attempting rescue from "
                  f"{len(leftover)} leftover segments "
                  f"(threshold: {args.rescue_min_length}bp)",
                  file=sys.stderr)

            rescued, rescued_groups, round_summaries = run_pass(
                segments, links, leftover, args, round_id=round_id)
            all_summaries.extend(round_summaries)

            n_rescued_comps = sum(1 for s in round_summaries
                                  if s["status"] == "ok")
            if not rescued:
                print(f"[round {round_id}] rescued 0 segments, stopping",
                      file=sys.stderr)
                break

            union_keep |= rescued
            node_group.update(rescued_groups)
            print(f"[round {round_id}] rescued {len(rescued)} segments "
                  f"across {n_rescued_comps} component(s)",
                  file=sys.stderr)
    else:
        if args.verbose:
            print("[rescue] disabled (--rescue-min-length 0)", file=sys.stderr)

    # ---- Write output. ----
    write_filtered_gfa(raw_lines, union_keep, node_group, args.gfa_out,
                       verbose=args.verbose)

    # ---- Tips file: TSV with one row per component (any round). ----
    # status: ok / too_few_tips / too_small / below_rescue_threshold
    tips_path = args.gfa_out + ".tips"
    with open(tips_path, "w") as tfh:
        tfh.write("\t".join([
            "comp_id", "round", "status", "size", "n_tips",
            "tip1", "tip2", "spine_len_bp",
            "ambiguous", "second_best_bp",
        ]) + "\n")
        for s in all_summaries:
            tfh.write("\t".join([
                str(s["comp_id"]),
                str(s["round"]),
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

    # ---- Summary stats. ----
    total_in = len(segments)
    total_kept = len(union_keep)
    n_ok = sum(1 for s in all_summaries if s["status"] == "ok")
    n_ambig = sum(1 for s in all_summaries if s["ambiguous"])
    n_rescued = sum(1 for s in all_summaries
                    if s["status"] == "ok" and s["round"] > 0)
    print(f"[done] {n_ok} component(s) had a spine identified "
          f"({n_ambig} ambiguous, {n_rescued} rescued)", file=sys.stderr)
    print(f"[done] kept {total_kept} / {total_in} segments "
          f"({100.0 * total_kept / max(total_in, 1):.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
