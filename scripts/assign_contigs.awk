#!/usr/bin/awk -f
#
# Assign assembly contigs to reference chromosomes from a MashMap PAF file.
#
# For each contig, sums alignment lengths per reference chromosome, then
# assigns the contig to the chromosome with the most aligned bases.
# Reports the runner-up chromosome and flags potential chimeras (where the
# second-best chromosome's aligned bases exceed CHIMERA_THRESHOLD * best).
#
# Usage:
#     ./assign_contigs.awk mashmap.paf > contig_assignments.tsv
#     awk -f assign_contigs.awk mashmap.paf > contig_assignments.tsv
#     awk -v CHIMERA_THRESHOLD=0.3 -f assign_contigs.awk mashmap.paf
#     awk -v MIN_FRAC=0.5 -f assign_contigs.awk mashmap.paf
#
# Options (set with -v):
#     CHIMERA_THRESHOLD   Flag contig as chimeric if 2nd-best chrom has at least
#                         this fraction of the best chrom's aligned bases. (0.2)
#     MIN_FRAC            Mark contigs as "low_coverage" if fraction of contig
#                         aligned to best chrom is below this. (0.3)
#
# Output columns (TSV):
#     contig             query contig name
#     contig_len         length of the contig (bp)
#     best_chrom         assigned reference chromosome
#     best_aln_bp        bases aligned to best_chrom
#     best_frac          best_aln_bp / contig_len
#     second_chrom       runner-up chromosome (or "." if none)
#     second_aln_bp      bases aligned to second_chrom (or 0)
#     second_ratio       second_aln_bp / best_aln_bp
#     flag               OK | CHIMERA | LOW_COVERAGE | CHIMERA,LOW_COVERAGE

BEGIN {
    FS = OFS = "\t"
    if (CHIMERA_THRESHOLD == "") CHIMERA_THRESHOLD = 0.2
    if (MIN_FRAC == "") MIN_FRAC = 0.3
}

# Sum alignment lengths per (contig, chromosome) pair.
# PAF columns: 1=query, 2=qlen, 3=qstart, 4=qend, 5=strand, 6=target, ...
{
    contig = $1
    qlen[contig] = $2
    chrom = $6
    aln_len = $4 - $3
    sum[contig SUBSEP chrom] += aln_len
    seen_chrom[contig SUBSEP chrom] = 1
    contigs[contig] = 1
}

END {
    # Header
    print "contig", "contig_len", "best_chrom", "best_aln_bp", "best_frac", \
          "second_chrom", "second_aln_bp", "second_ratio", "flag"

    for (c in contigs) {
        best_chrom = "unplaced"
        best_len = 0
        second_chrom = "."
        second_len = 0

        # Find best and second-best chromosome for this contig
        for (k in seen_chrom) {
            split(k, parts, SUBSEP)
            if (parts[1] != c) continue
            v = sum[k]
            if (v > best_len) {
                second_chrom = best_chrom
                second_len = best_len
                best_chrom = parts[2]
                best_len = v
            } else if (v > second_len) {
                second_chrom = parts[2]
                second_len = v
            }
        }

        # Compute fractions and flags
        contig_len = qlen[c]
        best_frac = (contig_len > 0) ? best_len / contig_len : 0
        second_ratio = (best_len > 0) ? second_len / best_len : 0

        flag = "OK"
        is_chimera = (second_len > 0 && second_ratio >= CHIMERA_THRESHOLD)
        is_lowcov = (best_frac < MIN_FRAC)
        if (is_chimera && is_lowcov) flag = "CHIMERA,LOW_COVERAGE"
        else if (is_chimera)         flag = "CHIMERA"
        else if (is_lowcov)          flag = "LOW_COVERAGE"

        printf "%s\t%d\t%s\t%d\t%.3f\t%s\t%d\t%.3f\t%s\n", \
            c, contig_len, best_chrom, best_len, best_frac, \
            second_chrom, second_len, second_ratio, flag
    }
}
