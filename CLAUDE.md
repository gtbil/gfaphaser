# Project notes for Claude Code

This file is loaded automatically by Claude Code on startup. It captures
context about the codebase that would otherwise need re-explaining each
session.

## What this crate does

Two binaries for processing hifiasm-style unitig GFA graphs:

- `gfa_split` — splits a GFA into per-connected-component files. Components
  whose unique sequence length is ≥ `--min-len` get their own file; smaller
  ones are concatenated into one `small_components.gfa`.
- `gfa_haps` — for each connected component, generates a pool of candidate
  walks via randomized bidirected traversal, then picks the best pair of
  walks (two haplotypes) under a lexicographic objective (maximize union
  coverage → maximize symmetric difference → minimize revisits). Runs a
  coverage-repair pass afterward to splice in any uncovered nodes.

See `README.md` for full usage.

## Project layout

```
Cargo.toml         # workspace manifest with one library and two binaries
src/lib.rs         # shared: parsing, bidirected Graph, DSU, Walk, Step
src/bin/gfa_split.rs
src/bin/gfa_haps.rs
```

The library exposes: `parse_gfa`, `Graph`, `Segment`, `Edge`, `Step`,
`Walk`, `DSU`, `cigar_overlap_len`, `connected_components`,
`entry_exit_for_orient`, `halfedge`. No external dependencies — pure
std-library Rust.

## Build & test

```
cargo build --release
./target/release/gfa_split <input.gfa> --min-len 100000
./target/release/gfa_haps <component.gfa> --candidates 128
```

Release build matters: candidate-pair scoring is O(K²) and benefits a lot
from optimization. Don't waste time on `cargo run` in debug mode for any
graph bigger than a few dozen nodes.

There are no unit tests yet — testing has been done by constructing
synthetic GFAs in Python and running both binaries end-to-end. If adding
tests, prefer integration tests that exercise the binaries on small
synthetic GFAs over unit tests of internal functions, because the
correctness criteria are mostly about emergent walk properties (coverage,
phasing, structural validity) rather than per-function behavior.

## Bidirected graph conventions (IMPORTANT — easy to get wrong)

Unitig graphs are bidirected. Each segment has two ends:

- `side 0` = the `-` end (left end of the unitig's forward sequence)
- `side 1` = the `+` end (right end)

An L-line `L A fo B to_o cigar` means: A's `fo`-side connects to B's
`to_o`-side, where:

- `fo == '+'` means A's right (side 1) attaches to the edge
- `fo == '-'` means A's left (side 0) attaches to the edge
- `to_o == '+'` means B's left (side 0) attaches — B is entered from the left
- `to_o == '-'` means B's right (side 1) attaches — B is entered from the right

So when parsing an L-line into adjacency:

- exit side of A = (1 if fo=='+' else 0)
- entry side of B = (0 if to_o=='+' else 1)

The reverse edge (traversing the same physical link from B back to A):
exit B on `b_entry_side`, enter A on `a_exit_side`. **The side does NOT
flip** — the edge attaches to specific physical ends. This was a real bug
I fixed early on; my first version had `1 - side` for the reverse, which
broke linear chains.

In a walk, each Step records an orientation char. The step's entry side
follows from orientation:

- step is '+' → entered on side 0, exits on side 1
- step is '-' → entered on side 1, exits on side 0

This is encoded in `entry_exit_for_orient()` in lib.rs.

## P-lines, not W-lines

We emit GFA P-lines (paths), not W-lines (walks). The GFA spec restricts
W-lines to overlap-free graphs; hifiasm unitig graphs have nonzero
overlaps on every L-line, so W-lines are spec-noncompliant. P-lines carry
an explicit Overlaps field listing the CIGAR for each consecutive
segment pair.

`Step::cigar_in` holds the original CIGAR string from the L-line that was
traversed to enter that step. It's preserved verbatim through parsing,
random walks, and the BFS used by coverage repair, so the P-line Overlaps
field is faithful to the input.

## Coverage repair architecture

After the best-pair selector picks two walks, `repair_coverage()` tries to
splice each uncovered node into one of the walks. Two splice flavors:

- `try_splice()`: BFS from walk step i's exit-state to target's exit-state,
  then BFS from target's exit-state to walk step j's entry-state. Replaces
  walk[i+1..=j] with the detour. **Must skip walk[j] when reattaching**
  because the detour ends at exactly walk[j]'s state — including it would
  produce adjacent duplicates. This was the bug behind the malformed walks
  Bandage rejected.
- `try_extend_end()`: BFS from walk's last exit-state to target, append.

Splices are only accepted if joint coverage **strictly increases** —
prevents flip-flop oscillation where splicing one node evicts another.

`is_valid_walk()` is run after every splice and again at output time.
It checks: (a) no adjacent same-node-same-orientation steps, (b) every
consecutive step pair is connected by a real L-line edge. If a malformed
walk ever sneaks through, it gets caught here and emits a loud error
rather than being written to the GFA.

## Things I considered but didn't do (yet)

- **Parallelism.** Per-component computation is embarrassingly parallel.
  Adding rayon would give a near-linear speedup on multi-component inputs.
- **Smarter complement walks.** The current "leader + complements that
  avoid leader's nodes" heuristic is OK but not optimal. A min-cost flow
  formulation could give provably optimal phasing on bubble chains.
- **Multi-node detour search.** `try_splice` handles one uncovered node at
  a time. For a cluster of unreachable-individually-but-reachable-together
  nodes, this can leave coverage on the table. A path-cover formulation
  would be cleaner.
- **A unit test harness.** Most "tests" are ad-hoc Python scripts that
  generate synthetic GFAs. Worth formalizing into `tests/` with cargo's
  integration testing.

## Style notes

- Prefer `eprintln!` for diagnostics, `println!` only for the
  tab-separated summary line consumed by downstream tooling.
- Don't add `clap` or other CLI crates unless really needed — hand-rolled
  arg parsing is intentional, keeps the crate dependency-free.
- I/O is line-oriented streaming. `gfa_haps` always loads sequences into
  memory (needed for FASTA output). `gfa_split` never needs them.
- The Rng (xorshift64*) is intentionally simple. If randomness quality
  ever matters, swap for the `rand` crate, but for now diversity comes
  from candidate count, not RNG quality.

## Recent history (for context)

In the conversation that produced this code:

1. Started with a single-file `gfa_split.rs`. Switched to a Cargo workspace
   with a shared lib when `gfa_haps` was added.
2. Found and fixed the bidirected-reverse-edge bug.
3. Added termination logic to walks so cycles don't make them loop forever
   on the revisit budget.
4. Added coverage repair with monotonic-coverage acceptance and an
   iteration cap.
5. Fixed an adjacent-duplicate bug in splice that produced Bandage-invalid
   walks; added `is_valid_walk` defense.
6. Switched W-lines to P-lines after user pointed out hifiasm graphs have
   overlaps and W-lines are spec-restricted.
