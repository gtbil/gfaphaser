# GFAphaser Pipeline — How-To Guide

## Overview

GFAphaser takes the primary unitig graph produced by **hifiasm** (`*.p_utg.gfa`) and separates it into two phased haplotype assemblies organized by chromosome. Hifiasm's unitig graphs capture both parental haplotypes in a single assembly graph — heterozygous regions appear as "bubbles" where two alternative paths represent the two alleles. This pipeline traverses that graph to extract one coherent path per haplotype, then maps the resulting sequences to a reference genome to produce chromosome-level FASTA files for each haplotype.

The final output is two FASTA files (`hap1.fa` and `hap2.fa`) with contigs named by chromosome, ready for downstream analysis.

### Algorithm overview

The pipeline operates in four stages. First, the unitig graph is normalized and split into connected components — one per chromosome or chromosome arm. Each component is then processed independently: alternating rounds of graph pruning remove structural tangles and spurious connections while preserving heterozygous bubble chains, leaving a clean backbone of the chromosome with its true heterozygous variation intact. The `gfa_haps` tool then walks this cleaned graph to find the pair of traversal paths that together maximizes the sequence covered, maximizes the sequence on which the two paths differ (driving them into opposite sides of heterozygous bubbles), and minimizes redundant revisits. A coverage repair pass splices in any nodes that neither walk visited. The two resulting paths — one per haplotype — are extracted as sequences, combined across all chromosomes, and finally mapped to a reference genome so that each contig can be labeled and organized by chromosome.

### Comparison with hifiasm's HiFi-only hap1/hap2 assemblies

Recent versions of hifiasm can output `hap1.p_ctg.gfa` and `hap2.p_ctg.gfa` from HiFi reads alone, without any additional data. However, without long-range phase information, hifiasm's haplotype assignment is **local**: each bubble in the graph is split into hap1 and hap2 sides independently, with no mechanism to stitch those local assignments into a globally consistent phase across the chromosome. The hap1 label at one bubble carries no guarantee of being the same parental chromosome as the hap1 label at the next bubble. The resulting hap1 and hap2 contig sets are therefore not truly complementary representations of two parental chromosomes — they are assemblies of locally-separated bubble halves that may switch parental identity at every intervening homozygous region.

GFAphaser addresses this by explicitly constructing two **end-to-end chromosome traversals**. Rather than assigning bubble halves independently, the algorithm walks the graph from one chromosome tip to the other and selects the pair of complete paths that maximizes the sequence on which they differ. Both haplotype walks are derived simultaneously, and the scoring objective drives them apart at every bubble they pass through together. The result is two sequences that are complementary by construction — wherever one walk takes the left side of a bubble, the other takes the right — assembled as coherent paths through the full chromosome rather than as a patchwork of independently labeled fragments.

**A note on switch errors:** a switch error is where a haplotype path inadvertently switches from representing one parental chromosome to the other partway through. GFAphaser is not immune to this: wherever a homozygous region (a node shared by both haplotypes) interrupts a bubble chain, both walks pass through the same unitig, and the algorithm has no signal to maintain phase consistency across that gap. The two haplotype identities can swap at each homozygous block. Hifiasm's HiFi-only assemblies face exactly the same limitation for exactly the same reason — the graph itself carries no long-range phase information across homozygous intervals. Switch errors should be understood as a fundamental property of sequence-graph-based phasing without external phase anchors, not as a defect specific to either tool.

---

## Prerequisites

### Input Files

| File | Location | Description |
|------|----------|-------------|
| `{ASM}.p_utg.gfa` | Parent directory of the working directory | Hifiasm primary unitig graph |
| `ref.fa` | Working directory | Reference genome for chromosome assignment |

### Required Software

All tools must be available on `$PATH`:

- `get_blunted` — cleans blunt-ended overlaps from the GFA
- `vg` — variation graph toolkit (graph normalization, format conversion, sequence extraction)
- `minimap2` — sequence alignment
- `mashmap` — fast reference mapping for chromosome assignment
- `miniasm` — overlap-based local assembly
- `samtools` — FASTA indexing and extraction

### Binaries from This Project

Build the Rust binaries first (release mode is required — scoring is slow in debug):

```bash
cargo build --release
cp target/release/gfa_split target/release/gfa_haps .
```

Ensure `prune_to_spine.py` and `prune_ports.py` (from `scripts/`) are accessible from the working directory.

### Cluster

The pipeline uses **SLURM** and submits to the `atlas` queue.

---

## Running the Pipeline

```bash
./scripts/submit_pipeline.sh <ASM> <PROJECT>
```

- `<ASM>` — your assembly identifier (e.g., `cotton_sample_1`); used to locate `{ASM}.p_utg.gfa` and name all outputs
- `<PROJECT>` — your SLURM account/project name for resource allocation

This submits all four stages as chained SLURM jobs. Stages 2–4 wait automatically for the prior stage to complete.

---

## Stage-by-Stage Walkthrough

### Stage 1 — Graph Normalization and Component Splitting

**Script:** `scripts/prep_gfa.sh` | **Resources:** 1 node, 48 cores, 374 GB RAM, 24 h max

**Removing overlaps with `get_blunted`**

In a hifiasm unitig graph, every link (L-line) records that the end of one unitig overlaps with the beginning of the next by some number of bases — these overlapping bases appear in the sequences of *both* adjacent unitigs. This is a faithful representation of how reads were assembled, but it creates a problem for any downstream graph operation: the same physical DNA bases are stored twice, once at the end of one node and again at the start of the next.

`get_blunted` converts the graph from this **overlap representation** to a **blunt representation**: it trims the overlapping bases from node ends so that each base in the genome appears in exactly one node. Edges become simple adjacency markers with no shared sequence. The `--provenance` flag writes a coordinate mapping file recording how original node positions translate to the new trimmed positions, which is needed to recover original coordinates later if necessary.

This conversion is a prerequisite for all subsequent graph operations. Without blunting:
- **Sequence extraction** would double-count the overlapping bases, producing incorrect contig lengths and sequences
- **Graph normalization** tools like `vg` expect blunt graphs — their compaction, unchop, and orientation algorithms assume each base lives in exactly one node
- **Path coordinates** would be ambiguous, since the same genomic position could be described as either the end of one node or the start of the next

After blunting, `vg` normalization compacts node IDs, breaks short cycles, orients all nodes forward, merges collinear nodes (unchop), and prunes tiny isolated subgraphs (< 50 kb). Normalization ensures a clean, consistently-oriented graph before splitting.

The normalized GFA is then split by **connected component** using `gfa_split`. Each connected component is a set of unitigs that are physically linked to each other — in a well-assembled diploid genome, each component typically corresponds to one chromosome or chromosome arm. Components with at least 100 kb of unique sequence each get their own file; smaller ones are combined into `small_components.gfa`.

**Outputs:**
```
subgraphs/{ASM}/gfa/{ASM}.component_01.gfa   # largest component
subgraphs/{ASM}/gfa/{ASM}.component_02.gfa
...
subgraphs/{ASM}/gfa/{ASM}.small_components.gfa
```

---

### Stage 2 — Haplotype Walk Extraction Per Component

**Script:** `scripts/make_haplotypes.sh` | **Resources:** 3 cores, 24 GB RAM, 12 h max per task  
**Runs as a SLURM array job** — one task per component, all running in parallel.

Each component is processed independently through a **9-step refinement pipeline** before haplotype extraction:

**Graph refinement (alternating passes):**

1. **`prune_to_spine.py`** — Identifies the chromosome backbone and removes everything not on it.

   A "tip" in the bidirected sense is a node where one port (left or right end) has zero incident edges — the sequence terminates there, with nothing to traverse further in that direction. These are the physical ends of chromosomes or chromosome arms.

   The algorithm works as follows:
   - Find all bidirected tips in the component.
   - Run Dijkstra's shortest-path algorithm from every tip to every other tip, weighting each step by the length (in bp) of the node being entered. This measures the base-pair distance between tips along the graph.
   - Select the tip pair with the greatest distance — the **diameter** of the component. This is the key insight for handling internal tangles: a dead-end node buried inside a tangle is graph-close to all surrounding tips because the tangle is locally dense; the true chromosome ends are graph-far apart because the full chromosome length separates them.
   - Find all nodes that lie on at least one simple path between the chosen tip pair. The method adds a virtual edge between the two chosen tips, then computes biconnected components of the resulting graph. Every node in the biconnected component containing that virtual edge lies on some simple path between the tips and is kept. Nodes outside it are discarded. Critically, **both sides of every heterozygous bubble between the tips are preserved** — each allele path is a valid simple path between the tips, so all bubble nodes pass this test.
   - A **rescue pass** then examines the discarded nodes. If the leftover graph contains any tip-to-tip component with a spine ≥ 100 kb, it is re-run through the same algorithm and saved as a separate component in the output. This handles secondary chromosome arms that were connected to the main chromosome through a tangle: once the tangle is removed, the arm becomes its own independent piece. Links between the original spine and any rescued components are severed.

2. **`prune_ports.py`** — Enforces that each node port has at most one incident edge, distinguishing legitimate bubbles from structural breaks.

   Each node in the blunted graph has two ports: a left port (L) and a right port (R). In a clean assembly graph that has been properly resolved by the assembler, each port should connect to at most one neighbor — the graph should consist of linear chains interrupted by bubbles, where each bubble diverges from one node, offers two alternative paths, and reconverges at another node. A port with more than one edge that does *not* reconverge is a structural artifact — typically a misassembly junction, a repeat-induced tangle, or a leftover connection from a different chromosome.

   The algorithm works as follows:
   - Identify every port with more than one incident edge.
   - For each such port, compute the **forward-reachable state set** from each edge leaving the port. Starting from a given (node, entered-port) state, a BFS follows the graph forward — exit via the opposite port, traverse all outgoing edges, enter each neighbor at their corresponding port — and collects every (node, entered-port) state reachable.
   - Two edges at the same port that share any state in their reachable sets will eventually converge on a common downstream node: this is a **legitimate bubble** representing a heterozygous region. All edges in such a group are kept.
   - Two edges whose reachable sets are completely disjoint will never converge — this is a **structural break**. The edges are clustered by a union-find over pairwise reconvergence, then all but the largest cluster are cut. Tiebreaking is by greedy chain length: the edge leading to the longest downstream path is preferred. Only edges (L-lines) are removed; all nodes (segments) are kept.

3. **VG normalization** between passes — re-applies compaction, orientation, and pruning to keep the graph clean after each structural change.

Three rounds of spine-pruning and port-pruning are applied, interspersed with normalization, progressively simplifying the graph until only the chromosome backbone and its bubble chains remain.

**Haplotype extraction (`gfa_haps`):**

Once the graph is clean, `gfa_haps` extracts two haplotype walks:

1. A pool of random **bidirected walks** is generated. A bidirected walk respects the physical orientation of each unitig — exiting a node on one end forces entry into the next node on a specific end, just as physical DNA strands do.

2. A second pool of walks is generated that is biased to diverge from the best single walk, encouraging the algorithm to cover alternative paths through bubbles (the two alleles).

3. The best **pair** of walks is selected using a lexicographic objective: first maximize the total sequence covered by either walk (union coverage), then maximize the sequence covered by exactly one walk (symmetric difference — this drives the two walks into opposite sides of bubbles), then minimize redundant revisits.

4. A **coverage repair** pass tries to splice any uncovered nodes back into one of the two walks, ensuring the output paths are as complete as possible.

Each component contributes two paths: `{ASM}_h1_component_NN` and `{ASM}_h2_component_NN`.

**Output per component:**
```
subgraphs/{ASM}/gfa_walks/{ASM}.{COMPONENT}.vg
```

---

### Stage 3 — Combine Components and Extract Sequences

**Script:** `scripts/finish_gfa.sh` | **Resources:** 1 node, 48 cores, 374 GB RAM, 24 h max

All per-component graph files are merged into a single combined graph. Paths shorter than 100 kb are filtered out (likely fragmented or artifactual components), and graph segments not covered by any surviving path are removed. The remaining sequences are extracted to FASTA and split by haplotype.

**Outputs in `subgraphs/{ASM}/gfa_walks/`:**
```
combined.vg        # merged graph
combined.gfa       # merged graph in GFA format
combined.fa        # all haplotype sequences
hap1.fa            # haplotype 1 sequences only
hap2.fa            # haplotype 2 sequences only
```

---

### Stage 4 — Chromosome Assignment and Local Assembly

**Script:** `scripts/split_by_ref.sh` | **Resources:** 1 node, 48 cores, 374 GB RAM, 24 h max

`hap1.fa` and `hap2.fa` are mapped to the reference genome using **MashMap** (minimum 50 kb segments, ≥ 90% identity). Each contig is assigned to its best-matching reference chromosome using `assign_contigs.awk`.

For each chromosome (A01–A13, D01–D13, PT, MT for the cotton subgenomes plus plastid/mitochondria), the pipeline:

1. Pulls out the contigs assigned to that chromosome for each haplotype
2. Uses **minimap2** to find overlaps between those contigs
3. Uses **miniasm** to perform a local overlap-based assembly, merging overlapping contigs into longer sequences
4. Appends any contigs that miniasm could not incorporate (non-overlapping sequences are kept as-is)
5. Renames all sequences with the chromosome prefix (e.g., `A01.h1.*`)

**Final outputs in `subgraphs/{ASM}/miniasm/`:**
```
hap1.fa            # all chromosomes, haplotype 1 — final assembly
hap2.fa            # all chromosomes, haplotype 2 — final assembly
hap1.fa.fai
hap2.fa.fai
hap1.assign.tsv    # contig → chromosome assignment table for hap1
hap2.assign.tsv    # contig → chromosome assignment table for hap2
A01/               # per-chromosome intermediate files
A02/
...
```

---

## Output Summary

| File | Description |
|------|-------------|
| `subgraphs/{ASM}/gfa/{ASM}.component_NN.gfa` | Normalized, split assembly graph components (Stage 1) |
| `subgraphs/{ASM}/gfa_walks/hap1.fa` | Full haplotype 1 assembly before chromosome splitting (Stage 3) |
| `subgraphs/{ASM}/gfa_walks/hap2.fa` | Full haplotype 2 assembly before chromosome splitting (Stage 3) |
| `subgraphs/{ASM}/miniasm/hap1.fa` | **Final** haplotype 1 assembly, organized by chromosome (Stage 4) |
| `subgraphs/{ASM}/miniasm/hap2.fa` | **Final** haplotype 2 assembly, organized by chromosome (Stage 4) |
| `subgraphs/{ASM}/miniasm/{CHR}/hap1.fa` | Chromosome-specific sequences for haplotype 1 |
| `subgraphs/{ASM}/miniasm/hap1.assign.tsv` | Contig-to-chromosome mapping table |
| `logs/out/` | SLURM stdout/stderr logs for all stages |

---

## Key Concepts

**GFA / unitig graph** — Graphical Fragment Assembly format. A hifiasm unitig graph represents the assembly as a set of nodes (unitigs — unambiguously assembled sequence segments) connected by edges (sequence overlaps). Heterozygous regions produce "bubbles" where two paths through the graph represent the two alleles.

**Connected component** — A subset of graph nodes where every node is reachable from every other node by following edges. In a diploid assembly, major chromosomes typically form separate connected components.

**Bidirected graph** — Each node in an assembly graph has two ends (the left end and the right end of the unitig sequence). Edges connect specific ends of nodes to specific ends of other nodes, capturing strand orientation. A valid walk through the graph must respect which end it enters and exits each node.

**Bubble** — A region of the graph where two paths diverge from a shared node and later reconverge. In a diploid organism, bubbles represent heterozygous variants — the two paths carry the two parental alleles. Preserving bubbles through graph refinement is essential for phasing.

**Coverage repair** — After selecting the best pair of haplotype walks, any graph nodes not visited by either walk are identified. The pipeline attempts to splice each uncovered node into one of the walks by finding a valid detour path, ensuring the output paths are as complete as possible.

**P-line** — GFA path record format used here (as opposed to W-lines). Hifiasm unitig graphs have nonzero overlaps on every edge, and the GFA specification restricts W-lines (walk records) to overlap-free graphs. P-lines include explicit overlap CIGAR strings for each consecutive segment pair and are used here for spec compliance.

---

## Troubleshooting and Tips

- **SLURM logs** are written to `logs/out/`. Check these first if a stage fails.
- **Low coverage on complex components:** Increase the walk pool size by editing `--candidates` in `make_haplotypes.sh` (default is 128 in the pipeline). Higher values improve coverage on graphs with many bubbles at the cost of runtime.
- **Component size threshold:** The `--min-len` argument to `gfa_split` (default 100 kb) controls which components get their own file vs. go into `small_components.gfa`. Adjust in `prep_gfa.sh` if you expect smaller chromosomes.
- **Performance:** Always use the release-mode binaries (`cargo build --release`). The candidate-pair scoring step is O(K²) in the number of candidate walks and is much too slow in debug mode for real-sized graphs.
- **Reference chromosomes:** The chromosome list in `split_by_ref.sh` (A01–A13, D01–D13, PT, MT) is hardcoded for cotton. Edit the `CHROMS` variable in that script if working with a different species.
