#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48                     # n processor core(s) per node X 2 threads per core
#SBATCH --mem=374GB               # maximum memory per node
#SBATCH -p atlas                  # standard node(s)
#SBATCH -J dotplot               # job name
#SBATCH --time 72:00:00
# Quick dotplot: assembly vs reference using nucmer + mummerplot
# Usage: ./dotplot.sh <reference.fasta> <assembly.fasta> [output_prefix]

source /home/${USER}/.bashrc


set -euo pipefail

module load mummer4

#REF="${1:?Usage: $0 <reference.fasta> <assembly.fasta> [prefix]}"
#QRY="${2:?Usage: $0 <reference.fasta> <assembly.fasta> [prefix]}"
#PREFIX=${REF}_vs_${QRY}
THREADS=48
#THREADS="${THREADS:-8}"

# 1. Align with nucmer (good defaults for whole-genome alignment)
#    --maxmatch  : use all anchor matches (best for repetitive genomes / divergent assemblies)
#    -l 100      : min exact match length (raise to 200+ for cleaner plots on close genomes)
#    -c 500      : min cluster length
nucmer -l 200 -c 1000 -b 500 -g 1000 -t "$THREADS" -p "$PREFIX" "$REF" "$QRY"

# 2. Filter alignments: keep 1-to-1 best mapping, min identity 90%, min length 1kb
delta-filter -1 -i 90 -l 5000 "${PREFIX}.delta" > "${PREFIX}.filter.delta"

# 3. Generate dotplot (PNG). --large for big genomes, --fat groups by sequence
mummerplot --png --large --fat \
    --filter --layout \
    -p "$PREFIX" \
    "${PREFIX}.filter.delta"

# 4. Also dump a coords table you can inspect / load into R/Python
show-coords -rclT "${PREFIX}.filter.delta" > "${PREFIX}.coords.tsv"

echo
echo "Done."
echo "  Dotplot : ${PREFIX}.png"
echo "  Coords  : ${PREFIX}.coords.tsv"
echo "  Delta   : ${PREFIX}.filter.delta"
