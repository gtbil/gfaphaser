#!/usr/bin/awk -f
#
# Filter a PAF file to keep only alignments where the target (reference) chrom
# matches the contig's assigned best_chrom from assign_contigs.awk output.
#
# Inputs (order matters!):
#   1. contig_assignments.tsv  (from assign_contigs.awk; pre-filter by flag yourself)
#   2. mashmap.paf             (the PAF to filter)
#
# Usage:
#     ./filter_paf_by_assignment.awk contig_assignments.tsv mashmap.paf > filtered.paf
#     awk -f filter_paf.awk contig_assignments.tsv mashmap.paf > filtered.paf
#
# Output: PAF lines from input #2 where column 1 (query/contig) is in the
# assignments table AND column 6 (target/chromosome) equals best_chrom.
#
# Assignments TSV columns expected (tab-separated, with header):
#   1=contig  2=contig_len  3=best_chrom  ... (rest ignored)

BEGIN { FS = OFS = "\t" }

# First file: build contig -> best_chrom map.
# Skip header (any line where the value of $2 is non-numeric).
NR == FNR {
    if ($2 ~ /^[0-9]+$/) {
        assign[$1] = $3
    }
    next
}

# Second file: emit PAF lines whose query is assigned AND target matches.
($1 in assign) && ($6 == assign[$1])
