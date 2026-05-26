#!/usr/bin/env python3
"""
paf_to_agp.py
=============

Build a reference-guided scaffolding AGP from a PAF of (query assembly) vs
(reference assembly) alignments.

For each query contig:
  1. Assign it to the reference chromosome where the most query bases align.
  2. Decide orientation by which strand carries more aligned query bases on
     that chromosome.
  3. Decide its position along the chromosome from the (length-weighted)
     median of target midpoints.
  4. Reject as ambiguous if the dominant chromosome accounts for less than
     --min-frac of all aligned query bases for the contig.

Then for each reference chromosome, sort placed contigs by position and emit
an AGP-1.1 record. Contigs are joined by N-gaps (default 100 bp). Unplaced
contigs are written as singleton scaffolds at the bottom under "chrUn".

Inputs
------
  --paf         PAF of query vs reference (FastGA, minimap2, etc).
                Query must be the *contigs you want to scaffold*; target must
                be the reference chromosomes. Query lengths are read from
                PAF column 2.
  --out         AGP output path (default stdout).
  --min-aln     Minimum alignment length to consider (bp). Default 5000.
  --min-frac    Minimum fraction of aligned bases on the winning chromosome
                for a contig to be placed (else -> unplaced). Default 0.70.
  --gap         Gap size between contigs in a scaffold (bp). Default 100.
  --report      Optional TSV summarizing the placement decision for each
                query contig.

  Note: contigs that don't appear in the PAF at all are skipped entirely.
        If you need to emit those as unplaced singletons too, you'd need to
        supply a separate contig list.

Output
------
  AGP-1.1 to stdout (or --out). The scaffolds are named after the reference
  chromosomes with a "_RagTag"-style suffix so downstream tools know these
  are reference-guided (we use "_scaff"); unplaced contigs go under "chrUn".
"""

import argparse
import sys
from collections import defaultdict
from statistics import median


def is_organelle(name):
    """True if the contig name looks like a mitochondrial or plastid sequence.

    Matches common naming conventions exactly (case-insensitive) rather than
    using substrings, to avoid false positives on names like 'A01' or
    contig IDs that happen to contain 'MT' or 'PT' as substrings.
    """
    n = name.upper()
    exact = {"MT", "PT", "M", "C", "CP",
             "CHRMT", "CHRPT", "CHRM", "CHRC", "CHRCP",
             "MITOCHONDRION", "MITOCHONDRIA", "PLASTID", "CHLOROPLAST"}
    if n in exact:
        return True
    # also match prefixes like "MT_*", "Pt_*", "chrMT_*", etc.
    prefixes = ("MT_", "PT_", "CHRMT_", "CHRPT_", "CHRM_", "CHRC_",
                "MITO_", "PLASTID_", "CHLORO_")
    return n.startswith(prefixes)


def parse_paf(path, min_aln):
    """Yield (qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend,
              n_matches, aln_len) for each line >= min_aln in *aln_len*.

    PAF columns (0-based):
      0 qname  1 qlen  2 qstart  3 qend  4 strand
      5 tname  6 tlen  7 tstart  8 tend
      9 nmatch 10 alnlen 11 mapq
    """
    with open(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2]); qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tlen = int(parts[6])
            tstart = int(parts[7]); tend = int(parts[8])
            aln_len = int(parts[10])
            if aln_len < min_aln:
                continue
            yield (qname, qlen, qstart, qend, strand,
                   tname, tlen, tstart, tend, aln_len)


def weighted_median(values, weights):
    """Median of *values* weighted by *weights* (both lists, same length)."""
    pairs = sorted(zip(values, weights))
    total = sum(w for _, w in pairs)
    cum = 0.0
    for v, w in pairs:
        cum += w
        if cum >= total / 2.0:
            return v
    return pairs[-1][0]


def detect_breakpoints(alns, window_bp, min_chunk_bp, min_chunk_frac):
    """Detect clean chimera breakpoints in a single query contig.

    Algorithm:
      1. Bin the contig into windows of `window_bp` along query coordinates.
      2. For each window, the "dominant" target is the chromosome with the
         most aligned query bases overlapping the window.
      3. Find runs of consecutive windows with the same dominant target.
      4. A run is a "chunk" if it spans at least `min_chunk_bp` and accounts
         for at least `min_chunk_frac` of the alignment bases in its span.
      5. If there are two or more chunks with different dominant targets,
         the contig is chimeric and is broken at the midpoints between
         adjacent chunks.

    `alns` is a list of dicts with keys: qstart, qend, strand, tname,
    tstart, tend, aln_len.

    Returns: a list of (qstart, qend, tname) tuples — the chunks to keep.
    If only one chunk survives, returns a single-element list (no split).
    Returns [] if no chunk is large enough (caller treats as unplaced).
    """
    if not alns:
        return []

    qlen = max(a["qend"] for a in alns)
    if qlen <= 0:
        return []

    # Build per-window dominant-target votes.
    n_windows = max(1, (qlen + window_bp - 1) // window_bp)
    # For each window, accumulate aligned bases per target.
    window_votes = [defaultdict(int) for _ in range(n_windows)]
    for a in alns:
        wstart = a["qstart"] // window_bp
        wend = min(n_windows - 1, (a["qend"] - 1) // window_bp)
        # distribute aligned bases across windows proportionally to overlap
        for w in range(wstart, wend + 1):
            w_lo = w * window_bp
            w_hi = (w + 1) * window_bp
            overlap = max(0, min(a["qend"], w_hi) - max(a["qstart"], w_lo))
            if overlap > 0:
                window_votes[w][a["tname"]] += overlap

    # Per-window dominant target (or None if window has no alignment)
    dominant = []
    for votes in window_votes:
        if not votes:
            dominant.append(None)
        else:
            dominant.append(max(votes.items(), key=lambda kv: kv[1])[0])

    # Find runs of identical dominant target. Treat None as "carry forward"
    # (don't break a chunk just because of an alignment gap).
    runs = []   # list of (start_window, end_window_inclusive, target)
    i = 0
    while i < n_windows:
        if dominant[i] is None:
            i += 1
            continue
        j = i
        while j + 1 < n_windows and (dominant[j + 1] is None or dominant[j + 1] == dominant[i]):
            j += 1
        runs.append((i, j, dominant[i]))
        i = j + 1

    # Convert runs to (qstart, qend, target) and filter by chunk size.
    chunks = []
    for (ws, we, tgt) in runs:
        cs = ws * window_bp
        ce = min(qlen, (we + 1) * window_bp)
        span = ce - cs
        if span < min_chunk_bp:
            continue
        # also require min_chunk_frac of the alignment bases in this span
        # to point at this target (filters out chunks with mixed signal)
        bases_total = 0
        bases_to_tgt = 0
        for a in alns:
            ov = max(0, min(a["qend"], ce) - max(a["qstart"], cs))
            if ov > 0:
                bases_total += ov
                if a["tname"] == tgt:
                    bases_to_tgt += ov
        if bases_total == 0:
            continue
        if bases_to_tgt / bases_total < min_chunk_frac:
            continue
        chunks.append((cs, ce, tgt))

    # Collapse adjacent surviving chunks that share the same target. This is
    # needed because a single anomalous window in the middle of a long
    # same-target stretch can split that stretch into two runs in the
    # `runs` pass above; once we filter the anomaly out for being too short,
    # the flanking same-target chunks shouldn't remain artificially split.
    merged = []
    for c in chunks:
        if merged and merged[-1][2] == c[2]:
            prev_cs, prev_ce, prev_tgt = merged[-1]
            merged[-1] = (prev_cs, c[1], prev_tgt)
        else:
            merged.append(c)
    chunks = merged

    # Also absorb "island" chunks: a chunk whose flankers BOTH map to the
    # same target, and which is itself small relative to its flankers, is
    # almost certainly a local repeat / duplication, not a real chimera.
    # Absorb it into one merged chunk of the flanker target. Repeat until
    # no more islands are found (a sweep may expose new opportunities).
    while True:
        new_merged = []
        changed = False
        i = 0
        while i < len(chunks):
            if (i + 2 < len(chunks)
                    and chunks[i][2] == chunks[i + 2][2]
                    and chunks[i + 1][2] != chunks[i][2]):
                left = chunks[i]
                island = chunks[i + 1]
                right = chunks[i + 2]
                left_size = left[1] - left[0]
                island_size = island[1] - island[0]
                right_size = right[1] - right[0]
                if island_size < min(left_size, right_size):
                    new_merged.append((left[0], right[1], left[2]))
                    i += 3
                    changed = True
                    continue
            new_merged.append(chunks[i])
            i += 1
        chunks = new_merged
        if not changed:
            break

    if len(chunks) <= 1:
        return chunks

    # If multiple chunks survive, adjust breakpoints to midpoints between
    # adjacent chunks rather than window boundaries, so the contig is
    # cleanly partitioned with no gaps and no overlaps.
    adjusted = []
    for k, (cs, ce, tgt) in enumerate(chunks):
        new_cs = 0 if k == 0 else (chunks[k - 1][1] + cs) // 2
        new_ce = qlen if k == len(chunks) - 1 else (ce + chunks[k + 1][0]) // 2
        adjusted.append((new_cs, new_ce, tgt))
    return adjusted


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out", default="-")
    ap.add_argument("--min-aln", type=int, default=5000)
    ap.add_argument("--min-frac", type=float, default=0.70)
    ap.add_argument("--gap", type=int, default=100)
    ap.add_argument("--report", default=None,
                    help="Optional TSV summarizing placement decisions")
    ap.add_argument("--split-chimeras", action="store_true",
                    help="Detect contigs aligning cleanly to multiple chromosomes "
                         "and break them into sub-contigs.")
    ap.add_argument("--split-window", type=int, default=500_000,
                    help="Window size (bp) for chimera detection. Default 500 kb.")
    ap.add_argument("--split-min-chunk", type=int, default=2_000_000,
                    help="Minimum chunk size (bp) to consider a contig segment as a "
                         "real placement, rather than noise. Default 2 Mb.")
    ap.add_argument("--split-min-frac", type=float, default=0.80,
                    help="Within a candidate chunk, fraction of aligned bases that "
                         "must point at the dominant chromosome. Default 0.80.")
    ap.add_argument("--overlap-min-bp", type=int, default=100_000,
                    help="Minimum target-coordinate overlap (bp) between two placed "
                         "contigs on the same chromosome to count as a real conflict. "
                         "Smaller overlaps are ignored as trivial. Default 100 kb.")
    ap.add_argument("--trim-min-keep", type=int, default=500_000,
                    help="After trimming an overlapping contig, the kept piece must "
                         "be at least this long (bp). Shorter remainders are dropped. "
                         "Default 500 kb.")
    ap.add_argument("--walks", default=None,
                    help="Optional path to write the per-chromosome reference "
                         "coverage table (same content shown in the log).")
    args = ap.parse_args()

    # --- First pass: collect raw alignments per contig and query lengths.
    # We keep raw alignments (not aggregated) so the chimera-detection step
    # can see per-alignment query coordinates. Aggregation by target happens
    # downstream, per unit (whole contig or chimera-split sub-contig).
    query_lengths = {}
    ref_lengths = {}                   # tname -> tlen, captured from PAF col 7
    contig_alns = defaultdict(list)   # qname -> list of alignment dicts

    for (qname, qlen, qstart, qend, strand,
         tname, tlen, tstart, tend, aln_len) in parse_paf(args.paf, args.min_aln):
        # PAF should report a consistent length per query; trust the first
        # and warn if a later record disagrees.
        prev = query_lengths.get(qname)
        if prev is None:
            query_lengths[qname] = qlen
        elif prev != qlen:
            print(f"WARNING: {qname} has inconsistent qlen in PAF "
                  f"({prev} vs {qlen}); keeping {prev}", file=sys.stderr)
        # Reference length: trust the first sighting per chromosome.
        ref_lengths.setdefault(tname, tlen)
        contig_alns[qname].append({
            "qstart": qstart, "qend": qend, "strand": strand,
            "tname": tname, "tstart": tstart, "tend": tend,
            "aln_len": aln_len,
        })

    # --- Decide placement per query contig (with optional chimera splitting) ---
    # placements is keyed by a unique "unit id". For un-split contigs the unit
    # id is the contig name; for split contigs it's "<contig>:start-end".
    # Each entry tracks the underlying contig and (start, end) slice.
    placements = {}   # unit_id -> dict
    unplaced = []     # list of unit_ids (unplaced contig names; we don't split
                      # unplaced contigs)
    chimera_split_count = 0

    def make_placement_from_alns(unit_id, source_contig, slice_start, slice_end, alns):
        """Compute a single placement decision from a list of alignments.
        The `alns` list is stored on the placement so that downstream
        trimming steps can re-evaluate the unit after coordinate changes."""
        targets = defaultdict(lambda: {
            "aln_bases": 0, "plus_bases": 0, "minus_bases": 0, "midpoints": [],
            "tstart_min": None, "tend_max": None,
        })
        for a in alns:
            t = targets[a["tname"]]
            t["aln_bases"] += a["aln_len"]
            if a["strand"] == "+":
                t["plus_bases"] += a["aln_len"]
            else:
                t["minus_bases"] += a["aln_len"]
            t["midpoints"].append(((a["tstart"] + a["tend"]) / 2.0, a["aln_len"]))
            if t["tstart_min"] is None or a["tstart"] < t["tstart_min"]:
                t["tstart_min"] = a["tstart"]
            if t["tend_max"] is None or a["tend"] > t["tend_max"]:
                t["tend_max"] = a["tend"]
        total = sum(t["aln_bases"] for t in targets.values())
        if not targets:
            return None
        best_chrom, best = max(targets.items(), key=lambda kv: kv[1]["aln_bases"])
        frac = best["aln_bases"] / total if total else 0.0
        strand = "+" if best["plus_bases"] >= best["minus_bases"] else "-"
        mids = [m for m, _ in best["midpoints"]]
        wts = [w for _, w in best["midpoints"]]
        pos = weighted_median(mids, wts)
        return {
            "unit_id": unit_id,
            "source_contig": source_contig,
            "slice_start": slice_start,   # 1-based, inclusive
            "slice_end": slice_end,       # 1-based, inclusive
            "chrom": best_chrom, "strand": strand, "pos": pos,
            "tstart": best["tstart_min"], "tend": best["tend_max"],
            "frac": frac, "total_aln": total, "status": "placed",
            "alns": list(alns),
        }

    for qname, alns in contig_alns.items():
        qlen = query_lengths[qname]

        # Optionally try to split this contig at chimeric breakpoints.
        chunks = None
        if args.split_chimeras:
            chunks = detect_breakpoints(
                alns, args.split_window, args.split_min_chunk, args.split_min_frac
            )

        # If splitting produced 2+ chunks, place each chunk as its own unit.
        if chunks is not None and len(chunks) >= 2:
            chimera_split_count += 1
            print(f"  chimera: splitting {qname} into {len(chunks)} pieces "
                  f"at q-coordinates {[(cs, ce, tgt) for cs, ce, tgt in chunks]}",
                  file=sys.stderr)
            for (cs, ce, _tgt_hint) in chunks:
                # gather just the alignments overlapping this slice
                sub_alns = [a for a in alns
                            if min(a["qend"], ce) - max(a["qstart"], cs) > 0]
                unit_id = f"{qname}:{cs + 1}-{ce}"
                p = make_placement_from_alns(unit_id, qname, cs + 1, ce, sub_alns)
                if p is None:
                    continue
                # apply min-frac threshold to each chunk independently
                if p["frac"] < args.min_frac or p["total_aln"] < args.min_aln:
                    p["status"] = "ambiguous"
                    placements[unit_id] = p
                    unplaced.append(unit_id)
                else:
                    placements[unit_id] = p
            continue

        # Otherwise, place the whole contig as a single unit.
        unit_id = qname
        p = make_placement_from_alns(unit_id, qname, 1, qlen, alns)
        if p is None:
            continue
        if p["frac"] < args.min_frac or p["total_aln"] < args.min_aln:
            p["status"] = "ambiguous"
            placements[unit_id] = p
            unplaced.append(unit_id)
        else:
            placements[unit_id] = p

    if args.split_chimeras:
        print(f"Chimeras detected and split: {chimera_split_count}", file=sys.stderr)

    # ------------------------------------------------------------------
    # Resolve target-coordinate overlaps between placed contigs.
    # Strategy:
    #   * If a contig B is fully enclosed in another (A.tstart <= B.tstart and
    #     B.tend <= A.tend on the same chromosome), drop B entirely — it has
    #     no unique reference territory.
    #   * Otherwise, trim B's query slice to remove the portion that aligned
    #     into A's target range. After trimming, re-compute B's placement
    #     from only the surviving alignments. If the trim splits B into two
    #     non-overlapping query pieces (A enclosed in B's target span), emit
    #     both as separate units.
    #   * When deciding "A trims B vs B trims A", the contig with more aligned
    #     bases on the shared chromosome is the trimmer (A); the other (B) is
    #     trimmed.
    #   * Iterate sweeps until no overlap above --overlap-min-bp remains, since
    #     trimming may shift target spans and create new conflicts.
    # ------------------------------------------------------------------

    def trim_unit_by_exclude(p, exclude_chrom, exclude_tstart, exclude_tend):
        """Return a list of new placement dicts produced by trimming `p` so
        that no surviving alignment on `exclude_chrom` falls inside
        [exclude_tstart, exclude_tend].

        Alignments to chromosomes other than `exclude_chrom` are not filtered
        by the exclude range (they live in a different coordinate namespace),
        but they only contribute to the new unit if they are spanned by the
        query slice of a surviving exclude-chromosome side.

        Algorithm:
          1. Partition `exclude_chrom` alignments into LEFT (tmid<exclude_tstart),
             RIGHT (tmid>exclude_tend), and DROP (tmid inside exclude).
          2. Each surviving side defines a query-coordinate bounding box.
          3. For each side, gather ALL alignments (any chrom) of `p` that fall
             inside that query bounding box, then re-run make_placement_from_alns.
          4. To guarantee progress, require that the new unit either places on
             a different chromosome OR its target span on `exclude_chrom`
             is strictly outside the exclude range. Otherwise drop the side.
        """
        on_chrom_left = []
        on_chrom_right = []
        for a in p["alns"]:
            if a["tname"] != exclude_chrom:
                continue
            tmid = (a["tstart"] + a["tend"]) / 2.0
            if tmid < exclude_tstart:
                on_chrom_left.append(a)
            elif tmid > exclude_tend:
                on_chrom_right.append(a)
            # else: midpoint inside exclude -> drop

        new_units = []
        sides = []
        if on_chrom_left:
            sides.append(("L", on_chrom_left))
        if on_chrom_right:
            sides.append(("R", on_chrom_right))

        for side_label, side_alns_on_chrom in sides:
            q_lo = min(a["qstart"] for a in side_alns_on_chrom)
            q_hi = max(a["qend"] for a in side_alns_on_chrom)
            q_lo = max(q_lo, p["slice_start"] - 1)
            q_hi = min(q_hi, p["slice_end"])
            if q_hi - q_lo < args.trim_min_keep:
                continue
            new_slice_start = q_lo + 1
            new_slice_end = q_hi
            # Re-gather ALL alignments (any chrom) whose query interval lies
            # entirely within the new slice. This is what defines the new
            # unit's placement context.
            sub_alns = [a for a in p["alns"]
                        if a["qstart"] >= q_lo and a["qend"] <= q_hi]
            if not sub_alns:
                continue
            uid = f"{p['source_contig']}:{new_slice_start}-{new_slice_end}"
            new_p = make_placement_from_alns(
                uid, p["source_contig"], new_slice_start, new_slice_end, sub_alns
            )
            if new_p is None:
                continue
            # Progress check: if the new unit still places on exclude_chrom AND
            # its target span still intersects the exclude range, drop this side.
            # Otherwise we'd loop forever.
            if (new_p["chrom"] == exclude_chrom
                    and new_p["tstart"] < exclude_tend
                    and new_p["tend"] > exclude_tstart):
                continue
            if new_p["frac"] < args.min_frac or new_p["total_aln"] < args.min_aln:
                new_p["status"] = "ambiguous"
            new_units.append(new_p)

        return new_units

    def find_one_conflict():
        """Find the first overlap conflict above --overlap-min-bp.
        Returns (trimmer_uid, trimmed_uid) or None if no conflict remains."""
        by_chr = defaultdict(list)
        for uid, p in placements.items():
            if p["status"] == "placed":
                by_chr[p["chrom"]].append(uid)
        for chrom, uids in by_chr.items():
            uids.sort(key=lambda u: placements[u]["tstart"])
            for i in range(len(uids)):
                pi = placements[uids[i]]
                for j in range(i + 1, len(uids)):
                    pj = placements[uids[j]]
                    if pj["tstart"] >= pi["tend"]:
                        break
                    overlap = min(pi["tend"], pj["tend"]) - max(pi["tstart"], pj["tstart"])
                    if overlap < args.overlap_min_bp:
                        continue
                    # winner = more aligned bases; loser is trimmed (or dropped
                    # entirely if enclosed)
                    if pi["total_aln"] >= pj["total_aln"]:
                        return uids[i], uids[j]
                    else:
                        return uids[j], uids[i]
        return None

    n_enclosed_dropped = 0
    n_trimmed = 0
    max_iters = 200    # safety net
    iters = 0
    while True:
        iters += 1
        if iters > max_iters:
            print(f"WARNING: overlap-resolution iteration cap ({max_iters}) reached", file=sys.stderr)
            break
        conflict = find_one_conflict()
        if conflict is None:
            break
        winner_uid, loser_uid = conflict
        a = placements[winner_uid]
        b = placements[loser_uid]
        # Enclosed check: loser b fully inside winner a -> drop b
        if a["tstart"] <= b["tstart"] and b["tend"] <= a["tend"]:
            print(f"  enclosed: {b['unit_id']} ({b['tstart']:,}-{b['tend']:,}) "
                  f"inside {a['unit_id']} ({a['tstart']:,}-{a['tend']:,}) on {a['chrom']} -> drop",
                  file=sys.stderr)
            b["status"] = "enclosed_dropped"
            if loser_uid not in unplaced:
                unplaced.append(loser_uid)
            n_enclosed_dropped += 1
            continue
        # Otherwise: trim loser to remove the overlapping segment
        new_units = trim_unit_by_exclude(b, a["chrom"], a["tstart"], a["tend"])
        # Mark the original loser as superseded (won't appear in AGP)
        b["status"] = "trimmed_out"
        if loser_uid not in unplaced:
            # only put it on the unplaced list if no trim survivors take its place
            pass
        # If trimming produced 0 surviving pieces, the loser becomes an
        # unplaced singleton.
        if not new_units:
            b["status"] = "trim_too_short"
            unplaced.append(loser_uid)
            print(f"  trim: {b['unit_id']} on {a['chrom']} fully consumed by overlap "
                  f"with {a['unit_id']}; no surviving piece long enough -> unplaced",
                  file=sys.stderr)
            continue
        n_trimmed += 1
        msg_pieces = []
        for nu in new_units:
            # avoid uid collisions: if a unit with this id already exists,
            # disambiguate by appending a small index. In practice the slice
            # coordinates make collisions extremely unlikely.
            uid = nu["unit_id"]
            suffix = 1
            while uid in placements:
                uid = f"{nu['unit_id']}.{suffix}"
                suffix += 1
            nu["unit_id"] = uid
            placements[uid] = nu
            if nu["status"] != "placed":
                unplaced.append(uid)
            msg_pieces.append(f"{uid} ({nu['tstart']:,}-{nu['tend']:,}, {nu['status']})")
        print(f"  trim: {b['unit_id']} -> {', '.join(msg_pieces)} "
              f"(overlap with {a['unit_id']} on {a['chrom']})",
              file=sys.stderr)

    if n_enclosed_dropped or n_trimmed:
        print(f"Overlap resolution: enclosed_dropped={n_enclosed_dropped}, "
              f"trimmed={n_trimmed}, iters={iters}", file=sys.stderr)

    # --- Group placed contigs by chromosome, sort by position ---
    by_chrom = defaultdict(list)
    for unit_id, p in placements.items():
        if p["status"] == "placed" and p["chrom"]:
            by_chrom[p["chrom"]].append((p["pos"], unit_id, p["strand"],
                                          p["source_contig"],
                                          p["slice_start"], p["slice_end"]))
    for chrom in by_chrom:
        by_chrom[chrom].sort()

    # --- Emit AGP ---
    out = sys.stdout if args.out == "-" else open(args.out, "w")
    print("##agp-version\t1.1", file=out)

    # one scaffold per reference chromosome; skip organelle chromosomes
    for chrom in sorted(by_chrom):
        if is_organelle(chrom):
            continue
        scaffold_name = chrom
        position = 1
        part = 0
        contigs_here = [c for c in by_chrom[chrom] if not is_organelle(c[3])]
        for i, (_, unit_id, strand, source_contig, sstart, send) in enumerate(contigs_here):
            seg_len = send - sstart + 1
            part += 1
            comp_start = position
            comp_end = position + seg_len - 1
            print("\t".join(map(str, [
                scaffold_name, comp_start, comp_end, part,
                "W", source_contig, sstart, send, strand
            ])), file=out)
            position = comp_end + 1

            if i < len(contigs_here) - 1:
                part += 1
                gap_start = position
                gap_end = position + args.gap - 1
                print("\t".join(map(str, [
                    scaffold_name, gap_start, gap_end, part,
                    "N", args.gap, "scaffold", "yes", "align_genus"
                ])), file=out)
                position = gap_end + 1

    # unplaced contigs -> singleton scaffolds named after the unit itself.
    # For chimera-split sub-contigs that ended up ambiguous, we still emit a
    # slice rather than the full contig. Skip units that were superseded by
    # trim products (status=trimmed_out) — they've been replaced.
    emitted_unplaced = set()
    for unit_id in sorted(unplaced):
        if unit_id in emitted_unplaced:
            continue
        p = placements[unit_id]
        if p["status"] == "trimmed_out":
            continue
        if is_organelle(p["source_contig"]):
            continue
        sstart = p["slice_start"]; send = p["slice_end"]
        seg_len = send - sstart + 1
        scaffold_name = unit_id
        print("\t".join(map(str, [
            scaffold_name, 1, seg_len, 1, "W", p["source_contig"], sstart, send, "+"
        ])), file=out)
        emitted_unplaced.add(unit_id)

    if args.out != "-":
        out.close()

    # --- Per-chromosome reference coverage summary ---
    # For each reference chromosome with at least one placed contig, report
    # the merged target-span coverage (non-overlapping union of placed
    # [tstart, tend) intervals) alongside the reference length.
    placed_count_by_chrom = defaultdict(int)
    intervals_by_chrom = defaultdict(list)
    for p in placements.values():
        if p["status"] != "placed" or p["chrom"] is None:
            continue
        if is_organelle(p["chrom"]):
            continue
        placed_count_by_chrom[p["chrom"]] += 1
        intervals_by_chrom[p["chrom"]].append((p["tstart"], p["tend"]))

    walked_by_chrom = {}
    for chrom, ivs in intervals_by_chrom.items():
        ivs.sort()
        merged_bp = 0
        cur_s, cur_e = ivs[0]
        for s, e in ivs[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged_bp += cur_e - cur_s
                cur_s, cur_e = s, e
        merged_bp += cur_e - cur_s
        walked_by_chrom[chrom] = merged_bp

    # Build rows once, then format for stderr and (optionally) write TSV.
    rows = []
    for chrom in sorted(walked_by_chrom):
        walked = walked_by_chrom[chrom]
        ref_len = ref_lengths.get(chrom, 0)
        pct = (100.0 * walked / ref_len) if ref_len else 0.0
        n_p = placed_count_by_chrom[chrom]
        ivs = intervals_by_chrom[chrom]
        first_bp = min(s for s, _ in ivs) + 1   # tstart is 0-based half-open
        last_bp = max(e for _, e in ivs)        # tend is 0-based half-open; last 1-based pos
        rows.append({
            "chrom": chrom, "placed": n_p,
            "first_bp": first_bp, "last_bp": last_bp,
            "walked_bp": walked, "ref_len": ref_len, "pct": pct,
        })

    total_walked = sum(r["walked_bp"] for r in rows)
    total_ref = sum(r["ref_len"] for r in rows)
    total_pct = (100.0 * total_walked / total_ref) if total_ref else 0.0

    # --- Formatted to stderr ---
    print("\nReference coverage by chromosome:", file=sys.stderr)
    print(f"{'chrom':<12}{'placed':>8}{'first_bp':>14}{'last_bp':>14}"
          f"{'walked_bp':>16}{'ref_len':>16}{'pct':>8}",
          file=sys.stderr)
    for r in rows:
        print(f"{r['chrom']:<12}{r['placed']:>8}{r['first_bp']:>14,}{r['last_bp']:>14,}"
              f"{r['walked_bp']:>16,}{r['ref_len']:>16,}{r['pct']:>7.1f}%",
              file=sys.stderr)
    if total_ref:
        print(f"{'TOTAL':<12}{'':>8}{'':>14}{'':>14}"
              f"{total_walked:>16,}{total_ref:>16,}{total_pct:>7.1f}%",
              file=sys.stderr)

    # --- TSV to --walks ---
    if args.walks:
        with open(args.walks, "w") as wf:
            wf.write("chrom\tplaced\tfirst_bp\tlast_bp\twalked_bp\tref_len\tpct\n")
            for r in rows:
                wf.write(f"{r['chrom']}\t{r['placed']}\t{r['first_bp']}\t"
                         f"{r['last_bp']}\t{r['walked_bp']}\t{r['ref_len']}\t"
                         f"{r['pct']:.4f}\n")
            if total_ref:
                wf.write(f"TOTAL\t\t\t\t{total_walked}\t{total_ref}\t"
                         f"{total_pct:.4f}\n")

    # --- Optional placement report ---
    if args.report:
        with open(args.report, "w") as rep:
            rep.write("unit\tsource_contig\tslice_start\tslice_end\t"
                      "status\tchrom\tstrand\tpos\tfrac\ttotal_aln_bp\n")
            for unit_id in sorted(placements):
                p = placements[unit_id]
                rep.write("\t".join([
                    unit_id,
                    p["source_contig"],
                    str(p["slice_start"]),
                    str(p["slice_end"]),
                    p["status"],
                    p["chrom"] or "",
                    p["strand"] or "",
                    f"{p['pos']:.0f}" if p["pos"] is not None else "",
                    f"{p['frac']:.4f}",
                    str(p["total_aln"]),
                ]) + "\n")

    # --- Summary to stderr ---
    n_placed = sum(1 for p in placements.values() if p["status"] == "placed")
    n_amb = sum(1 for p in placements.values() if p["status"] == "ambiguous")
    n_enc = sum(1 for p in placements.values() if p["status"] == "enclosed_dropped")
    n_tro = sum(1 for p in placements.values() if p["status"] == "trimmed_out")
    n_tts = sum(1 for p in placements.values() if p["status"] == "trim_too_short")
    print(f"Placed:    {n_placed}", file=sys.stderr)
    print(f"Ambiguous: {n_amb}", file=sys.stderr)
    if n_enc: print(f"Enclosed (dropped): {n_enc}", file=sys.stderr)
    if n_tro: print(f"Superseded by trim: {n_tro}", file=sys.stderr)
    if n_tts: print(f"Trim too short: {n_tts}", file=sys.stderr)


if __name__ == "__main__":
    main()
