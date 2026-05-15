# gfa_tools

Two binaries for slicing GFA unitig graphs into manageable pieces and
extracting two haplotype-style walks per piece. Pure Rust, no external
dependencies, builds with `cargo build --release`.

## Build

```
cargo build --release
# binaries land in ./target/release/{gfa_split, gfa_haps}
```

## `gfa_split` — split a GFA into per-component files

```
gfa_split <input.gfa> [--out-prefix PREFIX] [--min-len BP]
```

- Decomposes the graph into connected components by undirected link
  membership (orientation ignored for connectivity).
- Components with unique sequence length ≥ `--min-len` (default 100,000)
  each get their own file `<prefix>.component_NN.gfa`, indexed by
  descending size so `component_01` is the biggest.
- Smaller components are concatenated into `<prefix>.small_components.gfa`.
- "Unique sequence length" = Σ segment lengths − Σ link overlap CIGAR
  lengths. Exact for tree-shaped components; over-subtracts by one
  overlap per cycle-closing link in cyclic components.
- Routes `S`, `L`, `C`, `P`, `W` records into the right output by segment
  membership. `H` records are replicated to every output. `A` records are
  dropped. Each output is a self-contained valid GFA.

## `gfa_haps` — extract two haplotype walks per component

```
gfa_haps <input.gfa> [--endpoints tip|any] [--candidates K]
                     [--seed S] [--sample NAME] [--keep-sequences]
                     [--repair-max-revisits N] [--no-repair] [--verbose]
```

For each connected component, generates a diverse pool of candidate walks
and selects the **best pair** under a lexicographic objective:

1. Maximize **union coverage** (bp covered by either walk).
2. Maximize **symmetric difference** (bp covered by exactly one walk) —
   this is what drives the two walks to take opposite sides of bubbles.
3. Minimize **total revisits** across the two walks.

After pair selection, a **coverage repair** pass tries to splice any
uncovered nodes into one of the walks via shortest bidirected detours.
Splices are only accepted if they strictly increase joint coverage,
preventing oscillation. Use `--no-repair` to disable, or
`--repair-max-revisits N` (default 4) to cap the allowed revisits added
per splice.

### Endpoints (`--endpoints`)

- `tip` (default): walks start and end at graph tips — sides of nodes with
  no outgoing edges. If a component has no tips, falls back to `any`.
- `any`: walks can start and end anywhere.

### Candidates (`--candidates`)

Number of candidate walks generated per starting policy (default 64).
The actual pool is ~2K because the algorithm also generates K complementary
walks that strongly avoid the leader candidate's nodes (haplotype-divergence
pressure). All `~(2K)²/2` pairs are scored — cheap at typical K values.
Increase for hairier components if you suspect the best pair was missed.

### Walks are bidirected-aware

A walk respects unitig orientation: leaving a node via its right end forces
the next node to be entered on the side joined by the link. Walks emit GFA
W-line format with `>name` / `<name` for `+` / `-` orientations.

### Outputs

- `<basename>.with_paths.gfa` — the input GFA with two `P` lines added per
  component:
  ```
  P <sample>_h1_component_NN <seg1>+,<seg2>+,...,<segN>+ <cigar1>,<cigar2>,...
  P <sample>_h2_component_NN ...
  ```
  P-lines (not W-lines) are used because the GFA v1.1 spec restricts
  W-lines to graphs without overlaps between segments. Hifiasm unitig
  graphs have nonzero overlaps on every link, so P-lines are the correct
  format. CIGAR strings in the Overlaps field are preserved verbatim from
  the original L-lines.
- `<basename>.haps.fa` (only if `--keep-sequences`) — FASTA with one
  record per (component, haplotype), built by concatenating unitig
  sequences with overlaps trimmed from the prefix of each non-first step.

### A note on coverage

In a bidirected unitig graph, walks must traverse nodes coherently —
entering one side and exiting the other. This means dead-end side branches
(unitigs with no return path through their other side) cannot be visited
by a walk that also reaches a distal tip. With only two walks per
component, such branches can be unreachable. The tool reports
`nodes_covered / nodes_total` per component so you can spot this.

## Recommended pipeline

```
gfa_split input.gfa --min-len 100000
for comp in input.component_*.gfa; do
    gfa_haps "$comp" --keep-sequences --candidates 128
done
```

Per-component summary lines from `gfa_haps` are written to stdout as TSV,
suitable for piping into downstream analysis.
