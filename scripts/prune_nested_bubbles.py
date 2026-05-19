#!/usr/bin/env python3
"""
prune_nested_bubbles.py
=======================

Prune nodes that cause nested bubbles in a pangenome GFA, given a tip table
that names the anchor nodes the two haplotype paths must pass through.

Algorithm (per "ok" component in the tip table):
  1. Build the undirected graph of segments and links from the GFA, restricted
     to the segments in this component.
  2. Walk the tip "spine" t1 -> t2. (Your tip table currently lists at most two
     tips per component, so we have one bubble region to handle per component.
     The code generalises: if more tips appear in the future, it handles each
     consecutive pair.)
  3. Extract the bubble subgraph between consecutive tips: nodes that lie on
     ANY simple path from t_a to t_b. Implemented as: nodes reachable from t_a
     intersected with nodes that can reach t_b, after removing all OTHER tips
     so paths can't escape through them.
  4. Find the two node-disjoint simple paths from t_a to t_b whose combined bp
     of segment sequence is maximal. This is the max-weight 2-node-disjoint-
     paths problem; on a DAG it is solved exactly by min-cost max-flow with
     node capacities. If the bubble subgraph is not a DAG (back-edges within
     the bubble), we DAG-ify by orienting edges along a BFS from t_a; this is
     a heuristic but is the right behaviour for pangenome bubbles in practice.
  5. The union of nodes on those two paths is the KEEP set for this region.
     All other nodes in the bubble subgraph are pruned. We classify each pruned
     "blob" (connected component of pruned nodes) as either <50bp (small nested
     bubble — clean removal) or >=50bp (a genuine third alternative that lost
     the bp contest — also pruned, but logged separately so you can audit).
  6. Write a new GFA omitting pruned segments and any links touching them, plus
     a TSV report of what was removed and why.

Usage:
    python prune_nested_bubbles.py <gfa_in> <gfa_out> \\
        [--tips PATH] [--report PATH] \\
        [--small-bubble-bp 50] [--components 0,1]

The tip table is inferred as <gfa_in>.tips (override with --tips).
The report is inferred as <gfa_out> with its extension replaced by
'.prune_report.tsv' (override with --report).

Requirements: Python 3.8+, networkx >= 2.6.
    pip install networkx
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

try:
    import networkx as nx
except ImportError:
    sys.exit("This script requires networkx. Install with: pip install networkx")


# ---------------------------------------------------------------------------
# GFA parsing
# ---------------------------------------------------------------------------

@dataclass
class Segment:
    name: str
    seq: str
    length: int          # actual sequence length (len(seq)) OR LN tag if seq is '*'
    raw_line: str        # original line, so we can preserve all tags on write

@dataclass
class Link:
    from_seg: str
    from_orient: str     # '+' or '-'
    to_seg: str
    to_orient: str
    overlap: str         # e.g. '0M' — kept verbatim
    raw_line: str

def parse_gfa(path: Path) -> Tuple[Dict[str, Segment], List[Link], List[str]]:
    """
    Return (segments_by_name, links, other_lines).
    other_lines holds H/P/W/etc lines we just want to pass through (we will
    re-emit P and W lines only if all their segments survive — see write_gfa).
    """
    segments: Dict[str, Segment] = {}
    links: List[Link] = []
    others: List[str] = []

    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            tag = line[0]
            if tag == "S":
                parts = line.rstrip("\n").split("\t")
                # S <name> <seq> [tags...]
                name, seq = parts[1], parts[2]
                if seq == "*":
                    # Length must come from LN:i: tag
                    length = None
                    for t in parts[3:]:
                        if t.startswith("LN:i:"):
                            length = int(t[5:])
                            break
                    if length is None:
                        raise ValueError(
                            f"Segment {name} has seq='*' but no LN:i: tag — "
                            "cannot determine length")
                else:
                    length = len(seq)
                segments[name] = Segment(name, seq, length, line)
            elif tag == "L":
                parts = line.rstrip("\n").split("\t")
                # L <from> <fo> <to> <to> <ov> [tags...]
                links.append(Link(parts[1], parts[2], parts[3], parts[4],
                                  parts[5] if len(parts) > 5 else "*", line))
            else:
                others.append(line)

    return segments, links, others


# ---------------------------------------------------------------------------
# Tip table parsing
# ---------------------------------------------------------------------------

@dataclass
class TipRow:
    comp_id: int
    status: str
    size: int
    n_tips: int
    tips: List[str]      # only the named tips (tip1, tip2, ... — currently up to 2)
    spine_len_bp: int
    ambiguous: str
    second_best_bp: int

def parse_tips(path: Path) -> List[TipRow]:
    rows: List[TipRow] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Identify which columns hold tip names. The example header is:
        # comp_id status size n_tips tip1 tip2 spine_len_bp ambiguous second_best_bp
        tip_cols = [i for i, h in enumerate(header) if h.startswith("tip")]
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            tips = [f[i] for i in tip_cols if f[i] not in (".", "")]
            rows.append(TipRow(
                comp_id=int(f[idx["comp_id"]]),
                status=f[idx["status"]],
                size=int(f[idx["size"]]),
                n_tips=int(f[idx["n_tips"]]),
                tips=tips,
                spine_len_bp=int(f[idx["spine_len_bp"]]),
                ambiguous=f[idx["ambiguous"]],
                second_best_bp=int(f[idx["second_best_bp"]]),
            ))
    return rows


# ---------------------------------------------------------------------------
# Component assignment
# ---------------------------------------------------------------------------

def build_undirected(segments: Dict[str, Segment], links: List[Link]) -> nx.Graph:
    """Undirected adjacency for component finding + bubble subgraph extraction."""
    g = nx.Graph()
    for name in segments:
        g.add_node(name)
    for L in links:
        # We only care about node-level adjacency for the bubble extraction.
        # Orientation is preserved in the GFA on write; we don't need it for
        # finding paths between tips.
        if L.from_seg in segments and L.to_seg in segments:
            g.add_edge(L.from_seg, L.to_seg)
    return g

def find_components(g: nx.Graph) -> List[Set[str]]:
    """Return connected components as sets of node names, ordered by size desc.
    This matches the convention of most pangenome tools (comp 0 is the largest)."""
    comps = [set(c) for c in nx.connected_components(g)]
    comps.sort(key=len, reverse=True)
    return comps


# ---------------------------------------------------------------------------
# Bubble-region extraction and pruning
# ---------------------------------------------------------------------------

def bubble_region(g_und: nx.Graph,
                  comp_nodes: Set[str],
                  t_a: str,
                  t_b: str,
                  other_tips: Set[str]) -> Set[str]:
    """
    Nodes on some simple path from t_a to t_b within this component, where
    paths cannot pass through any tip other than t_a / t_b.

    We delete the other tips from the graph, find the nodes reachable from
    BOTH t_a and t_b in the remaining graph (i.e. on at least one path between
    them), and return that set including t_a and t_b.
    """
    sub_nodes = comp_nodes - other_tips
    if t_a not in sub_nodes or t_b not in sub_nodes:
        return set()
    h = g_und.subgraph(sub_nodes)
    reach_from_a = set(nx.node_connected_component(h, t_a))
    if t_b not in reach_from_a:
        return set()
    # All nodes in this connected piece of h that contain both t_a and t_b
    # are candidates for lying on some t_a–t_b path. (In a connected component
    # of an undirected graph, every node lies on some simple path between any
    # two given vertices iff removing it does not disconnect them — i.e. it is
    # not in a "side branch". For pruning we WANT to include side branches
    # rooted between t_a and t_b so we can prune them. The connected component
    # containing both tips is the right set.)
    return reach_from_a


def orient_as_dag(g_und: nx.Graph,
                  nodes: Set[str],
                  source: str) -> nx.DiGraph:
    """
    BFS from `source`; direct each edge from lower-BFS-depth endpoint to
    higher. Edges within the same BFS layer are dropped (they would create
    same-layer "siblings" with no defined direction — in a clean bubble they
    don't exist; if they do, they belong to a tangle we'd report).

    This is the standard DAG-ification used for bubble decomposition. It is a
    heuristic in the presence of true cycles, but pangenome bubble regions are
    near-acyclic by construction (the tips are chosen to be on both haplotypes),
    so any back-edges we encounter signal data weirdness, not a real cycle to
    preserve.
    """
    depth: Dict[str, int] = {source: 0}
    q = deque([source])
    while q:
        u = q.popleft()
        for v in g_und.neighbors(u):
            if v not in nodes:
                continue
            if v not in depth:
                depth[v] = depth[u] + 1
                q.append(v)

    dag = nx.DiGraph()
    dag.add_nodes_from(nodes)
    for u, v in g_und.subgraph(nodes).edges():
        du, dv = depth.get(u), depth.get(v)
        if du is None or dv is None:
            continue
        if du < dv:
            dag.add_edge(u, v)
        elif dv < du:
            dag.add_edge(v, u)
        # same depth → drop (would-be intra-layer edge)
    return dag


def two_disjoint_max_bp_paths(dag: nx.DiGraph,
                              source: str,
                              sink: str,
                              seg_len: Dict[str, int]
                              ) -> Optional[Tuple[List[str], List[str]]]:
    """
    Find the two node-disjoint s–t paths in `dag` that together maximise the
    sum of segment lengths of nodes used (each node counted at most once).

    Reduction to min-cost flow:
      - Split each internal node v into v_in -> v_out with capacity 1 and
        cost -seg_len[v] (negative because we want max bp).
      - Source and sink get capacity 2 (they are shared by both paths).
      - Each DAG edge u -> v becomes u_out -> v_in with capacity 1, cost 0.
      - Demand 2 units of flow from source_in to sink_out.

    Returns (path_a, path_b) as node lists, or None if two disjoint paths
    don't exist.
    """
    if source == sink or source not in dag or sink not in dag:
        return None

    flow = nx.DiGraph()
    IN, OUT = "_in", "_out"

    for v in dag.nodes():
        cap = 2 if v in (source, sink) else 1
        cost = -seg_len.get(v, 0)
        flow.add_edge(v + IN, v + OUT, capacity=cap, weight=cost)

    for u, v in dag.edges():
        flow.add_edge(u + OUT, v + IN, capacity=1, weight=0)

    # Super-source / super-sink supplying exactly 2 units.
    flow.add_node("SS", demand=-2)
    flow.add_node("ST", demand=2)
    flow.add_edge("SS", source + IN, capacity=2, weight=0)
    flow.add_edge(sink + OUT, "ST", capacity=2, weight=0)

    try:
        flow_dict = nx.min_cost_flow(flow)
    except nx.NetworkXUnfeasible:
        return None

    # Reconstruct two paths by following unit flow from source twice, deleting
    # used edges between extractions.
    def extract_one() -> Optional[List[str]]:
        path = [source]
        cur = source + OUT
        while not cur.startswith(sink):
            # Find an outgoing edge with positive flow
            nxt = None
            for w, units in flow_dict.get(cur, {}).items():
                if units > 0:
                    nxt = w
                    flow_dict[cur][w] -= 1
                    break
            if nxt is None:
                return None
            if nxt.endswith(IN):
                name = nxt[:-len(IN)]
                # consume the IN -> OUT split
                flow_dict[nxt][name + OUT] -= 1
                path.append(name)
                cur = name + OUT
            else:
                cur = nxt
        return path

    # Drain SS -> source_in edges twice
    flow_dict["SS"][source + IN] -= 2
    p1 = extract_one()
    if p1 is None:
        return None
    p2 = extract_one()
    if p2 is None:
        return None
    return p1, p2


# ---------------------------------------------------------------------------
# Pruning driver
# ---------------------------------------------------------------------------

@dataclass
class PruneRecord:
    comp_id: int
    tip_a: str
    tip_b: str
    bubble_nodes: int
    bubble_bp: int
    keep_nodes: int
    keep_bp: int
    pruned_blobs: List[Tuple[int, int, str]] = field(default_factory=list)
    # (n_nodes_in_blob, bp_in_blob, "small"|"large")
    warning: str = ""

def prune_component(comp_id: int,
                    comp_nodes: Set[str],
                    tips: List[str],
                    g_und: nx.Graph,
                    segments: Dict[str, Segment],
                    small_bubble_bp: int) -> Tuple[Set[str], List[PruneRecord]]:
    """
    Returns (pruned_node_set, list_of_records).
    """
    pruned: Set[str] = set()
    records: List[PruneRecord] = []
    seg_len = {n: segments[n].length for n in comp_nodes}

    for i in range(len(tips) - 1):
        t_a, t_b = tips[i], tips[i + 1]
        other_tips = set(tips) - {t_a, t_b}
        region = bubble_region(g_und, comp_nodes, t_a, t_b, other_tips)
        if not region or len(region) <= 2:
            # Nothing between the tips, or they aren't actually connected
            # within this bubble subgraph.
            continue

        bubble_bp_total = sum(seg_len[n] for n in region)

        dag = orient_as_dag(g_und, region, t_a)
        rec = PruneRecord(
            comp_id=comp_id, tip_a=t_a, tip_b=t_b,
            bubble_nodes=len(region), bubble_bp=bubble_bp_total,
            keep_nodes=0, keep_bp=0,
        )

        # Sanity: confirm t_b reachable in the DAG.
        if t_b not in nx.descendants(dag, t_a):
            rec.warning = ("tip_b not reachable from tip_a after DAG-ification; "
                           "graph may have an unusual structure here — skipping")
            records.append(rec)
            continue

        result = two_disjoint_max_bp_paths(dag, t_a, t_b, seg_len)
        if result is None:
            rec.warning = "could not find two node-disjoint paths between tips; skipping"
            records.append(rec)
            continue
        p1, p2 = result
        keep = set(p1) | set(p2)
        rec.keep_nodes = len(keep)
        rec.keep_bp = sum(seg_len[n] for n in keep)

        # Off-walk nodes in this bubble region get pruned. Group them into
        # connected blobs (in the undirected bubble subgraph minus the keep
        # set) for reporting.
        off = region - keep
        if off:
            off_sub = g_und.subgraph(off)
            for blob in nx.connected_components(off_sub):
                blob_bp = sum(seg_len[n] for n in blob)
                kind = "small" if blob_bp < small_bubble_bp else "large"
                rec.pruned_blobs.append((len(blob), blob_bp, kind))
            pruned |= off

        records.append(rec)

    return pruned, records


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_gfa(out_path: Path,
              segments: Dict[str, Segment],
              links: List[Link],
              others: List[str],
              pruned: Set[str]) -> None:
    """
    Write a GFA omitting pruned segments and any L/P/W lines that reference them.
    H lines and other non-S/L are passed through unchanged unless they're P or W.
    """
    with out_path.open("w") as out:
        for line in others:
            tag = line[0]
            if tag in ("P", "W"):
                # P: <name> <segs comma-list> <overlaps>
                # W: <sample> <hap> <seq> <start> <end> <walk>
                # If any referenced segment was pruned, drop the line (we don't
                # try to patch paths — the user can re-derive them).
                parts = line.rstrip("\n").split("\t")
                if tag == "P":
                    seg_field = parts[2]
                    seg_names = [s.rstrip("+-") for s in seg_field.split(",")]
                else:  # W
                    walk_field = parts[6] if len(parts) > 6 else ""
                    # Walks look like ">1<2>3"; split on > and <
                    seg_names = []
                    cur = ""
                    for ch in walk_field:
                        if ch in "><":
                            if cur:
                                seg_names.append(cur)
                            cur = ""
                        else:
                            cur += ch
                    if cur:
                        seg_names.append(cur)
                if any(s in pruned for s in seg_names):
                    continue  # drop
                out.write(line)
            else:
                out.write(line)

        for name, seg in segments.items():
            if name in pruned:
                continue
            out.write(seg.raw_line if seg.raw_line.endswith("\n")
                      else seg.raw_line + "\n")

        for L in links:
            if L.from_seg in pruned or L.to_seg in pruned:
                continue
            out.write(L.raw_line if L.raw_line.endswith("\n")
                      else L.raw_line + "\n")


def write_report(report_path: Path, records: List[PruneRecord]) -> None:
    with report_path.open("w") as out:
        out.write("\t".join([
            "comp_id", "tip_a", "tip_b",
            "bubble_nodes", "bubble_bp",
            "kept_nodes", "kept_bp",
            "pruned_blob_count",
            "pruned_small_blobs", "pruned_small_bp",
            "pruned_large_blobs", "pruned_large_bp",
            "warning",
        ]) + "\n")
        for r in records:
            small = [b for b in r.pruned_blobs if b[2] == "small"]
            large = [b for b in r.pruned_blobs if b[2] == "large"]
            out.write("\t".join(str(x) for x in [
                r.comp_id, r.tip_a, r.tip_b,
                r.bubble_nodes, r.bubble_bp,
                r.keep_nodes, r.keep_bp,
                len(r.pruned_blobs),
                len(small), sum(b[1] for b in small),
                len(large), sum(b[1] for b in large),
                r.warning,
            ]) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa_in", type=Path, help="Input GFA file")
    ap.add_argument("gfa_out", type=Path, help="Output pruned GFA")
    ap.add_argument("--tips", type=Path, default=None,
                    help="Tip table TSV. Default: <gfa_in>.tips")
    ap.add_argument("--report", type=Path, default=None,
                    help="Output TSV report. Default: <gfa_out> with extension "
                         "replaced by '.prune_report.tsv'")
    ap.add_argument("--small-bubble-bp", type=int, default=50,
                    help="Threshold for classifying a pruned blob as 'small' "
                         "(reporting only — does not change what is pruned). "
                         "Default 50.")
    ap.add_argument("--components", type=str, default=None,
                    help="Comma-separated component IDs to process. "
                         "Default: every row with status='ok' in the tip table.")
    args = ap.parse_args()

    tips_path = args.tips or args.gfa_in.with_suffix(args.gfa_in.suffix + ".tips")
    report_path = args.report or args.gfa_out.with_suffix(".prune_report.tsv")

    print(f"[1/5] Reading GFA: {args.gfa_in}", file=sys.stderr)
    segments, links, others = parse_gfa(args.gfa_in)
    print(f"      {len(segments)} segments, {len(links)} links", file=sys.stderr)

    print(f"[2/5] Reading tips: {tips_path}", file=sys.stderr)
    tip_rows = parse_tips(tips_path)

    if args.components:
        wanted = set(int(c) for c in args.components.split(","))
        tip_rows = [r for r in tip_rows if r.comp_id in wanted]
    else:
        tip_rows = [r for r in tip_rows if r.status == "ok"]
    print(f"      processing {len(tip_rows)} component(s): "
          f"{[r.comp_id for r in tip_rows]}", file=sys.stderr)

    print("[3/5] Building component index", file=sys.stderr)
    g_und = build_undirected(segments, links)
    comps = find_components(g_und)
    # Map comp_id (ordinal by descending size, as in the tip table) to node set
    comp_lookup = {i: nodes for i, nodes in enumerate(comps)}

    print("[4/5] Pruning bubble regions", file=sys.stderr)
    all_pruned: Set[str] = set()
    all_records: List[PruneRecord] = []
    for row in tip_rows:
        if row.comp_id not in comp_lookup:
            print(f"  WARN comp {row.comp_id} not found in GFA — skipping",
                  file=sys.stderr)
            continue
        comp_nodes = comp_lookup[row.comp_id]
        # Verify tips are actually in this component; if not, the tip table
        # and GFA disagree about component numbering. Bail loudly.
        missing = [t for t in row.tips if t not in comp_nodes]
        if missing:
            print(f"  WARN comp {row.comp_id}: tips {missing} not in component "
                  f"(size {len(comp_nodes)}). Component ordering between GFA "
                  "and tip table may differ. Skipping.", file=sys.stderr)
            continue
        pruned, recs = prune_component(
            row.comp_id, comp_nodes, row.tips, g_und,
            segments, args.small_bubble_bp,
        )
        all_pruned |= pruned
        all_records.extend(recs)
        print(f"  comp {row.comp_id}: pruned {len(pruned)} nodes "
              f"({sum(segments[n].length for n in pruned)} bp)", file=sys.stderr)

    print(f"[5/5] Writing outputs", file=sys.stderr)
    write_gfa(args.gfa_out, segments, links, others, all_pruned)
    write_report(report_path, all_records)

    total_bp_pruned = sum(segments[n].length for n in all_pruned)
    print(f"\nDone. Pruned {len(all_pruned)} segments totalling {total_bp_pruned} bp.",
          file=sys.stderr)
    print(f"Output GFA: {args.gfa_out}", file=sys.stderr)
    print(f"Report:     {report_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
