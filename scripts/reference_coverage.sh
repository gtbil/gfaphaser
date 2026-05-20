#!/usr/bin/env bash
#
# reference_coverage.sh
#
# Identify which segments of the reference genome are covered by assembly
# contigs (via a MashMap PAF) and which are not. Useful for spotting candidate
# gap regions before running a gap-filler.
#
# Outputs (in OUTDIR):
#     ref_covered.bed       merged reference intervals covered by contigs
#     ref_uncovered.bed     reference intervals NOT covered by any contig
#     ref_uncovered.filt.bed  uncovered intervals >= MIN_GAP_SIZE bp
#     coverage_summary.tsv  per-chromosome covered/uncovered totals & %
#
# Usage:
#     ./reference_coverage.sh -p mashmap.paf -f reference.fasta [-o outdir] [-m 1000]
#
# Requires: bedtools, samtools (for .fai)

set -euo pipefail

# Defaults
OUTDIR="ref_coverage_out"
MIN_GAP_SIZE=1000

usage() {
    cat <<EOF
Usage: $0 -p <mashmap.paf> -f <reference.fasta> [-o <outdir>] [-m <min_gap_bp>]

Required:
  -p PAF        MashMap PAF output (assembly contigs mapped to reference)
  -f FASTA      Reference genome FASTA (needs .fai or will be created)

Optional:
  -o OUTDIR     Output directory (default: $OUTDIR)
  -m MIN_GAP    Minimum gap size in bp to report (default: $MIN_GAP_SIZE)
  -h            Show this help

Outputs:
  ref_covered.bed       merged covered reference intervals
  ref_uncovered.bed     all uncovered reference intervals
  ref_uncovered.filt.bed  uncovered intervals >= MIN_GAP bp
  coverage_summary.tsv  per-chromosome covered/uncovered bp and %

Example:
  $0 -p assembly_vs_ref.paf -f reference.fasta -o gaps_out -m 5000
EOF
    exit 1
}

PAF=""
REF=""

while getopts "p:f:o:m:h" opt; do
    case $opt in
        p) PAF="$OPTARG" ;;
        f) REF="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        m) MIN_GAP_SIZE="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "$PAF" || -z "$REF" ]] && usage
[[ ! -f "$PAF" ]] && { echo "Error: PAF file not found: $PAF" >&2; exit 1; }
[[ ! -f "$REF" ]] && { echo "Error: Reference FASTA not found: $REF" >&2; exit 1; }

# Check dependencies
for tool in bedtools samtools awk sort; do
    command -v $tool >/dev/null 2>&1 || { echo "Error: $tool not found in PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

# --- Step 1: ensure reference .fai exists, then make chromsizes file ---
if [[ ! -f "${REF}.fai" ]]; then
    echo "Indexing reference..." >&2
    samtools faidx "$REF"
fi

CHROMSIZES="$OUTDIR/ref.chromsizes"
cut -f1,2 "${REF}.fai" > "$CHROMSIZES"

# --- Step 2: PAF -> reference-coordinate BED ---
# PAF columns: 1=qname 2=qlen 3=qstart 4=qend 5=strand 6=tname 7=tlen 8=tstart 9=tend
# We want BED of reference (target) intervals: tname, tstart, tend
echo "Converting PAF to reference BED..." >&2
PAF_BED="$OUTDIR/ref_alignments.bed"
awk 'BEGIN{OFS="\t"} {print $6, $8, $9, $1, $4-$3, $5}' "$PAF" \
    | sort -k1,1 -k2,2n > "$PAF_BED"

# --- Step 3: merge overlapping reference intervals -> covered ---
echo "Merging overlapping covered intervals..." >&2
COVERED="$OUTDIR/ref_covered.bed"
bedtools merge -i "$PAF_BED" > "$COVERED"

# --- Step 4: complement -> uncovered intervals ---
echo "Computing uncovered regions..." >&2
UNCOVERED="$OUTDIR/ref_uncovered.bed"
# bedtools complement needs sorted chromsizes matching the BED's chrom order
sort -k1,1 "$CHROMSIZES" > "$OUTDIR/ref.chromsizes.sorted"
sort -k1,1 -k2,2n "$COVERED" > "$OUTDIR/ref_covered.sorted.bed"
bedtools complement -i "$OUTDIR/ref_covered.sorted.bed" -g "$OUTDIR/ref.chromsizes.sorted" \
    > "$UNCOVERED"

# --- Step 5: filter uncovered by minimum size ---
UNCOVERED_FILT="$OUTDIR/ref_uncovered.filt.bed"
awk -v m="$MIN_GAP_SIZE" 'BEGIN{OFS="\t"} ($3-$2) >= m {print $1,$2,$3,$3-$2}' \
    "$UNCOVERED" > "$UNCOVERED_FILT"

# --- Step 5b: coverage multiplicity (depth) across reference ---
# bedtools genomecov -bga gives BedGraph including 0-coverage intervals.
# Columns: chrom, start, end, depth
echo "Computing coverage multiplicity..." >&2
sort -k1,1 "$CHROMSIZES" > "$OUTDIR/ref.chromsizes.sorted"
MULT_BG="$OUTDIR/ref_multiplicity.bedgraph"
bedtools genomecov -bga -i "$PAF_BED" -g "$OUTDIR/ref.chromsizes.sorted" > "$MULT_BG"

# Intervals with depth >= 2 are "over-covered" -- candidate duplications,
# repeat collapses in the ref, or unresolved haplotypes in your assembly.
MULT_HIGH="$OUTDIR/ref_multiplicity.ge2.bed"
awk -v m="$MIN_GAP_SIZE" 'BEGIN{OFS="\t"}
    $4 >= 2 && ($3-$2) >= m {print $1,$2,$3,$3-$2,$4}' "$MULT_BG" > "$MULT_HIGH"

# --- Step 6: per-chromosome summary ---
echo "Building summary..." >&2
SUMMARY="$OUTDIR/coverage_summary.tsv"

# Build summary in a single awk pass; portable (no gawk extensions).
# We pass four files and dispatch by FILENAME.
awk -v fchrom="$CHROMSIZES" -v fcov="$COVERED" -v fgap="$UNCOVERED_FILT" -v fmult="$MULT_BG" '
     BEGIN{
        OFS="\t"
        print "chrom","chrom_len","covered_bp","uncovered_bp","pct_covered",
              "n_gaps_filt","largest_gap","mean_depth","max_depth","bp_depth_ge2","pct_depth_ge2"
     }
     FILENAME == fchrom { clen[$1] = $2; order[++n] = $1; next }
     FILENAME == fcov   { cov[$1] += $3 - $2; next }
     FILENAME == fgap   {
         ngaps[$1]++
         sz = $3 - $2
         if (sz > maxg[$1]) maxg[$1] = sz
         next
     }
     FILENAME == fmult  {
         w = $3 - $2
         d = $4
         depth_sum[$1] += d * w     # for mean depth (length-weighted)
         if (d > maxd[$1]) maxd[$1] = d
         if (d >= 2) bp_ge2[$1] += w
         next
     }
     END {
         for (i = 1; i <= n; i++) {
             c = order[i]
             L = clen[c]
             cb = (c in cov) ? cov[c] : 0
             ub = L - cb
             pct = (L > 0) ? 100.0 * cb / L : 0
             ng = (c in ngaps) ? ngaps[c] : 0
             mg = (c in maxg) ? maxg[c] : 0
             md = (c in depth_sum && L > 0) ? depth_sum[c] / L : 0
             mx = (c in maxd) ? maxd[c] : 0
             bp2 = (c in bp_ge2) ? bp_ge2[c] : 0
             pct2 = (L > 0) ? 100.0 * bp2 / L : 0
             printf "%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%.3f\t%d\t%d\t%.2f\n", \
                 c, L, cb, ub, pct, ng, mg, md, mx, bp2, pct2
         }
     }' "$CHROMSIZES" "$COVERED" "$UNCOVERED_FILT" "$MULT_BG" > "$SUMMARY"

# --- Tidy ---
rm -f "$OUTDIR/ref.chromsizes.sorted" "$OUTDIR/ref_covered.sorted.bed"

# --- Report ---
echo "" >&2
echo "=== Done ===" >&2
echo "Output directory: $OUTDIR" >&2
echo "" >&2
echo "Files:" >&2
echo "  $COVERED" >&2
echo "  $UNCOVERED" >&2
echo "  $UNCOVERED_FILT  (gaps >= ${MIN_GAP_SIZE} bp)" >&2
echo "  $MULT_BG         (per-base coverage depth, BedGraph)" >&2
echo "  $MULT_HIGH       (depth >= 2 intervals >= ${MIN_GAP_SIZE} bp)" >&2
echo "  $SUMMARY" >&2
echo "" >&2

# Quick stats to stderr
total_ref=$(awk '{s+=$2} END{print s}' "$CHROMSIZES")
total_cov=$(awk '{s+=$3-$2} END{print s+0}' "$COVERED")
total_unc=$(awk '{s+=$3-$2} END{print s+0}' "$UNCOVERED")
n_gaps=$(wc -l < "$UNCOVERED_FILT")
bp_ge2=$(awk '$4>=2 {s+=$3-$2} END{print s+0}' "$MULT_BG")
max_depth=$(awk 'BEGIN{m=0} $4>m{m=$4} END{print m}' "$MULT_BG")
echo "Reference total:     $(printf "%'d" $total_ref) bp" >&2
echo "Covered (>=1x):      $(printf "%'d" $total_cov) bp ($(awk -v a=$total_cov -v b=$total_ref 'BEGIN{printf "%.2f", 100*a/b}')%)" >&2
echo "Uncovered:           $(printf "%'d" $total_unc) bp" >&2
echo "Depth >= 2x:         $(printf "%'d" $bp_ge2) bp ($(awk -v a=$bp_ge2 -v b=$total_ref 'BEGIN{printf "%.2f", 100*a/b}')%)" >&2
echo "Max depth observed:  ${max_depth}x" >&2
echo "Gaps >= ${MIN_GAP_SIZE} bp:     $n_gaps" >&2
