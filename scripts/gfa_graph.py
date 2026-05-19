"""
gfa_graph.py
------------
Shared GFA loader + bounded forward-walk used by both
gfa_to_kmer_fasta.py and annotate_gfa_with_kmc_depth.py.

A "node end" is the pair (segment_id, orientation, side) where side is 'R'
(the right / 3' end as emitted in that orientation) or 'L' (the left / 5'
end). Outgoing links from the right end of (S, +) are the GFA L lines with
A=S, oA=+; outgoing links from the right end of (S, -) are those with A=S,
oA=-; etc.  In a GFA, a link L A oA B oB means "after traversing A in
direction oA you can next traverse B in direction oB."

walk_forward(start_seg, start_orient, target_len, ...) yields all paths of
total length up to target_len bases beginning at the right end of
(start_seg, start_orient). Each yielded path is a string of bases of length
min(target_len, total_available_path_length).  Output is truncated to
exactly target_len when longer.

Cycles are naturally bounded by target_len (each step adds >=1 base), so
recursion always terminates.

Branching can explode in tangled regions; max_paths caps total emitted
paths per call (default 4096) and yields a sentinel-free truncated set,
setting a flag on the supplied list `overflow_flag` if hit.
"""

from collections import defaultdict


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def oriented_seq(seq, orient):
    """Return seq as emitted when the segment is traversed in `orient`."""
    return seq if orient == "+" else revcomp(seq)


def flip(orient):
    return "+" if orient == "-" else "-"


def load_gfa(path):
    """
    Parse a GFA. Returns (seqs, out_links) where:
        seqs[node_id] = uppercased sequence string
        out_links[(node_id, orient)] = list of (next_node_id, next_orient)

    A link `L A oA B oB ...` means after (A,oA) you can go to (B,oB).
    By GFA semantics, this also implies after (B, flip(oB)) you can go to
    (A, flip(oA)) -- traversing the reverse complement path. We add both.
    """
    seqs = {}
    out_links = defaultdict(list)

    with open(path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t", 3)
                nid, seq = parts[1], parts[2]
                if seq != "*":
                    seqs[nid] = seq.upper()
            elif line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                a, oa, b, ob = parts[1:5]
                out_links[(a, oa)].append((b, ob))
                out_links[(b, flip(ob))].append((a, flip(oa)))

    # Dedup links to avoid double-emission when a GFA already includes both
    # directions of a link explicitly.
    for key in list(out_links.keys()):
        out_links[key] = list(set(out_links[key]))

    return seqs, out_links


def walk_forward(start_seg, start_orient, target_len, seqs, out_links,
                 max_paths=4096, overflow=None):
    """
    Yield distinct path-strings of length up to target_len bases that
    extend forward from the right end of (start_seg, start_orient).

    Each yielded string is built by concatenating successor segments'
    oriented sequences. If a path can be extended to >= target_len bases,
    it is truncated to exactly target_len.

    If more than max_paths distinct paths exist, only the first max_paths
    are yielded and overflow[0] is set True (if overflow is a 1-element list).
    """
    yielded = set()
    n_emitted = 0

    # Iterative DFS using an explicit stack of (current_node, current_orient,
    # accumulated_string).
    stack = [(start_seg, start_orient, "")]

    while stack:
        node, orient, acc = stack.pop()
        successors = out_links.get((node, orient), [])

        if not successors:
            # Dead end -- emit whatever we accumulated (could be empty).
            out = acc
            if out and out not in yielded:
                yielded.add(out)
                yield out
                n_emitted += 1
                if n_emitted >= max_paths:
                    if overflow is not None:
                        overflow[0] = True
                    return
            continue

        for nxt_node, nxt_orient in successors:
            if nxt_node not in seqs:
                continue
            nxt_bases = oriented_seq(seqs[nxt_node], nxt_orient)
            new_acc = acc + nxt_bases

            if len(new_acc) >= target_len:
                out = new_acc[:target_len]
                if out not in yielded:
                    yielded.add(out)
                    yield out
                    n_emitted += 1
                    if n_emitted >= max_paths:
                        if overflow is not None:
                            overflow[0] = True
                        return
            else:
                # Keep extending.
                stack.append((nxt_node, nxt_orient, new_acc))
