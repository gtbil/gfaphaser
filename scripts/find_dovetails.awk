#!/usr/bin/awk -f
#
# Find dovetail-style end-to-end overlaps between contigs from a minimap2 PAF.
# A dovetail is an alignment that touches one end of the query AND one end
# of the target, suggesting the two contigs can be joined with a small trim.
#
# Input: PAF from minimap2 -x asm5 -X (all-vs-all contigs, no self hits).
#
# Usage:
#     ./find_dovetails.awk contigs_ava.paf > dovetails.tsv
#
# Options (set with -v):
#     END_SLACK     bp of slack allowed from contig end to call it "at end" (500)
#     MIN_OVERLAP   minimum overlap length in bp (1000)
#     MIN_IDENT     minimum percent identity 0-100 (95)
#     DEDUP         if 1, only report one orientation per contig pair (1)
#
# PAF columns used:
#   1=qname 2=qlen 3=qstart 4=qend 5=strand 6=tname 7=tlen 8=tstart 9=tend
#   10=matches 11=aln_len 12=mapq
#
# Output columns (TSV):
#   q_contig  q_len  q_end_used  t_contig  t_len  t_end_used  strand
#   overlap_len  pct_identity  q_trim  t_trim  join_type  q_aln_range  t_aln_range
#
# end_used values: 5p, 3p
# join_type: how the two contigs would be joined after trimming
#   5p-5p, 5p-3p, 3p-5p, 3p-3p
# q_trim: bp to trim off the q end being used (negative means extend, positive means trim)
# t_trim: bp to trim off the t end being used

BEGIN {
    FS = OFS = "\t"
    if (END_SLACK == "")   END_SLACK = 500
    if (MIN_OVERLAP == "") MIN_OVERLAP = 1000
    if (MIN_IDENT == "")   MIN_IDENT = 95
    if (DEDUP == "")       DEDUP = 1

    print "q_contig","q_len","q_end_used","t_contig","t_len","t_end_used","strand",
          "overlap_len","pct_identity","q_trim","t_trim","join_type",
          "q_aln_range","t_aln_range"
}

# Skip self-hits in case -X wasn't used
$1 == $6 { next }

{
    qname = $1; qlen = $2; qstart = $3; qend = $4
    strand = $5
    tname = $6; tlen = $7; tstart = $8; tend = $9
    matches = $10; aln_len = $11

    # Identity: matches / aln_len (PAF doesn't carry exact identity, this is the
    # standard approximation. Some minimap2 versions emit cg/de tags too.)
    pct_id = (aln_len > 0) ? 100.0 * matches / aln_len : 0
    if (pct_id < MIN_IDENT) next

    # Overlap length: use mean of query and target alignment spans
    q_span = qend - qstart
    t_span = tend - tstart
    overlap = (q_span + t_span) / 2
    if (overlap < MIN_OVERLAP) next

    # --- Determine which end of the query the alignment hits ---
    # 5' = start (qstart near 0); 3' = end (qend near qlen)
    q_at_5p = (qstart <= END_SLACK)
    q_at_3p = (qlen - qend <= END_SLACK)

    # --- Same for target ---
    t_at_5p = (tstart <= END_SLACK)
    t_at_3p = (tlen - tend <= END_SLACK)

    # Must hit at least one end on each contig
    if (!(q_at_5p || q_at_3p)) next
    if (!(t_at_5p || t_at_3p)) next

    # Pick the end actually used. If both ends are within slack (very small
    # contig fully contained in alignment), prefer the closer end.
    if (q_at_5p && q_at_3p)      q_end_used = (qstart <= (qlen - qend)) ? "5p" : "3p"
    else if (q_at_5p)            q_end_used = "5p"
    else                         q_end_used = "3p"

    if (t_at_5p && t_at_3p)      t_end_used = (tstart <= (tlen - tend)) ? "5p" : "3p"
    else if (t_at_5p)            t_end_used = "5p"
    else                         t_end_used = "3p"

    # --- Reject containments ---
    # If both ends of query are inside the alignment, query is contained in target.
    # If both ends of target are inside the alignment, target is contained in query.
    # These aren't dovetails; they suggest one contig is a fragment of the other.
    if (q_at_5p && q_at_3p && (qend - qstart) >= qlen - 2*END_SLACK) {
        # query fully covered -> containment
        next
    }
    if (t_at_5p && t_at_3p && (tend - tstart) >= tlen - 2*END_SLACK) {
        next
    }

    # --- Compute trim amounts ---
    # q_trim: extra bp on the q end past the alignment (positive = trim this much)
    if (q_end_used == "5p") q_trim = qstart
    else                    q_trim = qlen - qend

    if (t_end_used == "5p") t_trim = tstart
    else                    t_trim = tlen - tend

    join_type = q_end_used "-" t_end_used

    # --- Sanity check on orientation ---
    # On + strand, a valid dovetail joins q-3p to t-5p (q_end overlaps t_start)
    #                     or q-5p to t-3p (q_start overlaps t_end)
    # On - strand, the target end interpretation flips because the target is
    # rev-comp'd relative to the query. In PAF, tstart/tend are still in target
    # forward coords, but the alignment direction is reversed.
    # We accept all combinations here and let the user inspect; for canonical
    # dovetails on +, you'd want 3p-5p or 5p-3p.

    # --- Deduplicate ---
    # The same overlap typically appears twice in all-vs-all PAF (A vs B and B vs A).
    # Keep only one by canonicalizing the contig pair name order.
    if (DEDUP) {
        if (qname < tname) key = qname "\t" tname
        else               key = tname "\t" qname
        if (key in seen) next
        seen[key] = 1
    }

    q_range = qstart "-" qend
    t_range = tstart "-" tend

    printf "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%.2f\t%d\t%d\t%s\t%s\t%s\n",
           qname, qlen, q_end_used,
           tname, tlen, t_end_used,
           strand,
           overlap, pct_id,
           q_trim, t_trim,
           join_type,
           q_range, t_range
}
