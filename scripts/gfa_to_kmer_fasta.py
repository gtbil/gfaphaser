#!/usr/bin/env python3
"""
gfa_to_kmer_fasta.py
--------------------
Emit a FASTA enumerating every k-mer POSITION in a GFA, handling segments
shorter than k correctly via bounded multi-segment forward walks.

For each segment S we emit:
    1. The full sequence of S in '+' orientation (handles all k-mer positions
       interior to S, and boundary-crossing positions on segments long enough
       to contain their own anchor).
    2. For each orientation oS in {+, -}, the string:
            anchor(S, oS) + fwd_path
       where anchor(S, oS) is the last (k-1) bases of S as emitted in oS
       (or all of S if shorter than k-1), and fwd_path is each distinct
       (k-1)-base extension reachable through outgoing links.
       These records contain exactly the k-mer positions starting in the
       last k-1 positions of (S, oS), which is precisely the set of
       boundary-crossing positions.

We emit walks in BOTH orientations because a k-mer position starting near
the END of (S, +) is a *different* graph position from one starting near
the END of (S, -). Both are real, and KMC will canonicalize them to the
same canonical k-mer when counting, so multiplicity comes out right.

Cycles are bounded by k-1 base accumulation, so termination is guaranteed.
Branching is capped per-anchor at MAX_PATHS to prevent combinatorial blowup;
overflow triggers a warning.

Usage:
    python gfa_to_kmer_fasta.py <input.gfa> <k> <output.fa>

Then build the graph KMC DB:
    kmc -k<k> -ci1 -cs1000000 -fm <output.fa> <graph_db_prefix> <tmpdir>
"""

import sys

from gfa_graph import (
    load_gfa,
    oriented_seq,
    walk_forward,
)


MAX_PATHS_PER_ANCHOR = 4096


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    gfa_in, k_str, fa_out = sys.argv[1:]
    k = int(k_str)
    overhang = k - 1

    seqs, out_links = load_gfa(gfa_in)
    print(f"[info] loaded {len(seqs)} segments, "
          f"{sum(len(v) for v in out_links.values())} link-ends",
          file=sys.stderr)

    n_segs_out = 0
    n_segs_short = 0
    n_walks_out = 0
    n_overflow_anchors = 0

    with open(fa_out, "w") as fout:
        # Pass 1: emit segments themselves. Skip segments shorter than k
        # since they hold zero interior k-mer windows; their boundary-crossing
        # k-mers are emitted by the Pass-2 walks from upstream segments.
        for nid, seq in seqs.items():
            if len(seq) < k:
                n_segs_short += 1
                continue
            fout.write(f">S_{nid}\n{seq}\n")
            n_segs_out += 1

        # Pass 2: emit boundary-crossing walks from each (segment, orient).
        for nid in seqs:
            for orient in ("+", "-"):
                anchor_full = oriented_seq(seqs[nid], orient)
                anchor = anchor_full[-overhang:] if len(anchor_full) >= overhang else anchor_full

                # If there are no outgoing links, there is nothing to emit
                # beyond what the segment itself already covered.
                if not out_links.get((nid, orient)):
                    continue

                overflow_flag = [False]
                n_for_this_anchor = 0
                for fwd in walk_forward(
                    nid, orient, overhang, seqs, out_links,
                    max_paths=MAX_PATHS_PER_ANCHOR,
                    overflow=overflow_flag,
                ):
                    if not fwd:
                        continue
                    rec = anchor + fwd
                    if len(rec) < k:
                        # No k-mer windows in this record; skip emission to
                        # keep the FASTA small.
                        continue
                    fout.write(f">W_{nid}{orient}_{n_for_this_anchor}\n{rec}\n")
                    n_walks_out += 1
                    n_for_this_anchor += 1

                if overflow_flag[0]:
                    n_overflow_anchors += 1
                    if n_overflow_anchors <= 5:
                        print(f"[warn] forward-walk path cap "
                              f"({MAX_PATHS_PER_ANCHOR}) hit at "
                              f"({nid},{orient}); graph multiplicity will be "
                              "approximate in this region.", file=sys.stderr)

    if n_overflow_anchors:
        print(f"[warn] {n_overflow_anchors} anchor(s) hit the path cap total.",
              file=sys.stderr)
    print(f"[done] wrote {n_segs_out} segments and {n_walks_out} walk records "
          f"(skipped {n_segs_short} segments shorter than k={k}) "
          f"-> {fa_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
