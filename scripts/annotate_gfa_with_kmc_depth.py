#!/usr/bin/env python3
"""
annotate_gfa_with_kmc_depth.py
------------------------------
Annotate every S line of a GFA with fair-share k-mer depth from a KMC
database, using a SORTED `odgi kmers` output for per-node attribution.

Memory is O(num_nodes), independent of kmers file size.

The kmers file is expected to be sorted lexicographically by canonical
k-mer (column 1 after canonicalization). The recommended pre-processing
pipeline is:

    odgi kmers -i graph.og -k 29 -c \\
        | awk 'BEGIN{OFS="\\t"} {
              # Emit canonical(kmer)\\tnode_id
              k=$1; n=split($2,a,":"); node=a[1]
              # Reverse-complement k
              rc=""; for(i=length(k);i>0;i--){
                  c=substr(k,i,1)
                  rc = rc ( c=="A"?"T":c=="T"?"A":c=="C"?"G":c=="G"?"C":"N" )
              }
              ck = (k<rc?k:rc)
              print ck, node
          }' \\
        | LC_ALL=C sort -k1,1 -S 8G --parallel 8 -T $TMPDIR \\
        > graph.kmers.sorted

The awk is a bit fiddly; an equivalent Python preprocessor is fine too.
What matters is the output is two tab-separated columns
(canonical_kmer, node_id) sorted by column 1.

Inputs:
    reads_db          KMC database (native format), prefix without extensions
    sorted_kmers      Sorted two-column TSV: canonical_kmer\\tnode_id
    gfa_in            Input GFA
    gfa_out           Annotated output GFA

For each canonical k-mer (one contiguous run of lines), we do ONE KMC
lookup, then distribute `read_count / run_length` to each node in the run.

Tags added to each S line:
    KC:i:<rounded_attributed_sum>   GFA standard k-mer count tag
    DP:f:<median_attributed_depth>  per-position median (robust)
    ad:f:<mean_attributed_depth>    per-position mean
    kn:i:<n_positions>              number of k-mer positions home to this node

Existing KC / DP / dp / ad / kn tags are stripped before writing.

Usage:
    python annotate_gfa_with_kmc_depth.py \\
        <reads_db_prefix> <sorted_kmers> <in.gfa> <out.gfa>
"""

import statistics
import sys
import time
from collections import defaultdict

import py_kmc_api as pka


STRIP_TAGS = {"KC", "DP", "dp", "ad", "kn"}


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    reads_prefix, sorted_path, gfa_in, gfa_out = sys.argv[1:]

    reads_db = pka.KMCFile()
    if not reads_db.OpenForRA(reads_prefix):
        sys.exit(f"ERROR: could not open reads KMC DB '{reads_prefix}'.")
    k = reads_db.KmerLength()
    print(f"[info] KMC k = {k}", file=sys.stderr)

    kmer_buf = pka.KmerAPI(k)
    rc_buf = pka.Count()

    # ---- Streaming attribution. ---------------------------------------------
    print("[info] streaming sorted kmers ...", file=sys.stderr)
    t0 = time.time()

    pools = defaultdict(list)

    cur_kmer = None
    cur_nodes = []
    n_runs = 0
    n_positions = 0
    n_bad = 0

    def flush():
        nonlocal n_runs
        if cur_kmer is None:
            return
        kmer_buf.from_string(cur_kmer)
        rc = rc_buf.value if reads_db.CheckKmer(kmer_buf, rc_buf) else 0
        mult = len(cur_nodes)
        contrib = rc / mult
        for nid in cur_nodes:
            pools[nid].append(contrib)
        n_runs += 1

    with open(sorted_path) as fh:
        for line in fh:
            try:
                ck, node_id = line.rstrip("\n").split("\t", 1)
            except ValueError:
                n_bad += 1
                continue
            if len(ck) != k:
                n_bad += 1
                continue
            n_positions += 1

            if ck != cur_kmer:
                flush()
                cur_kmer = ck
                cur_nodes = [node_id]
            else:
                cur_nodes.append(node_id)
        flush()

    t1 = time.time()
    if n_bad:
        print(f"[warn] {n_bad} malformed lines skipped", file=sys.stderr)
    print(f"[info] {n_positions} k-mer positions, {n_runs} unique canonical "
          f"k-mers, {len(pools)} nodes attributed in {t1 - t0:.1f}s",
          file=sys.stderr)

    # ---- Rewrite the GFA. ---------------------------------------------------
    t0 = time.time()
    n_segments = 0
    n_missing = 0
    with open(gfa_in) as fin, open(gfa_out, "w") as fout:
        for line in fin:
            if not line.startswith("S\t"):
                fout.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            node_id, seq = parts[1], parts[2]
            existing_tags = [t for t in parts[3:]
                             if t.split(":", 1)[0] not in STRIP_TAGS]

            pool = pools.get(node_id)
            if not pool:
                if len(seq) >= k:
                    n_missing += 1
                new_tags = ["KC:i:0", "DP:f:0", "ad:f:0", "kn:i:0"]
            else:
                total = sum(pool)
                median = statistics.median(pool)
                mean = total / len(pool)
                new_tags = [
                    f"KC:i:{int(round(total))}",
                    f"DP:f:{median:.4f}",
                    f"ad:f:{mean:.4f}",
                    f"kn:i:{len(pool)}",
                ]

            fout.write("\t".join(["S", node_id, seq] + existing_tags + new_tags)
                       + "\n")
            n_segments += 1

    t1 = time.time()
    if n_missing:
        print(f"[warn] {n_missing} segments of length >= k had no k-mers in "
              "the sorted kmers file. Check that the GFA and kmers file are "
              "from the same graph.", file=sys.stderr)
    print(f"[done] annotated {n_segments} segments in {t1 - t0:.1f}s "
          f"-> {gfa_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
