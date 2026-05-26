#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=300G
#SBATCH -t 4:00:00
#SBATCH -J fastga_aln
#SBATCH -A gbru_malvaceae
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

# Submit with:
#   GENOMES=mydir OUTPUT_PREFIX=run1 sbatch --export=ALL --array=1-$N fastga_aln.sh
# Defaults: GENOMES=genomes, OUTPUT_PREFIX=. (current directory)
#
# Setup before submitting:
#   GENOMES=${GENOMES:-genomes}
#   OUTPUT_PREFIX=${OUTPUT_PREFIX:-.}
#   mkdir -p $OUTPUT_PREFIX/{logs,alns,chained,plots,beds,pafs}
#   ls $GENOMES/ | grep '\.1gdb$' | grep -v '^Jin668v1\.1gdb$' > query_list.txt
#   N=$(wc -l < query_list.txt)
#   GENOMES=$GENOMES OUTPUT_PREFIX=$OUTPUT_PREFIX \
#       sbatch --export=ALL \
#              -o "$OUTPUT_PREFIX/logs/fastga_aln_%A_%a.out" \
#              -e "$OUTPUT_PREFIX/logs/fastga_aln_%A_%a.err" \
#              --array=1-$N fastga_aln.sh
#
# Note: Slurm directives are read at submit time and don't see env vars, so
# -o/-e are passed on the sbatch command line above to put logs under
# $OUTPUT_PREFIX/logs/. If you skip that override, logs land in cwd.

set -euo pipefail

# Inputs and outputs are parameterized via env vars.
GENOMES="${GENOMES:-genomes}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-.}"

# Reference (target) is fixed; query is read from the list (bare filenames).
ref="${GENOMES}/Jin668v1.1gdb"
qry="${GENOMES}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" query_list.txt)"

na=$(basename "${ref%.1gdb}")
nb=$(basename "${qry%.1gdb}")

echo "[task ${SLURM_ARRAY_TASK_ID}] aligning ${nb} (query) vs ${na} (target)"
echo "  GENOMES=${GENOMES}  OUTPUT_PREFIX=${OUTPUT_PREFIX}"

FastGA -T32 -i0.90 -l5000 -f5 -P"${TMPDIR}" \
    -1:${OUTPUT_PREFIX}/alns/${nb}_vs_${na}.1aln "$qry" "$ref"

ALNchain -c0.80 -e0.001 \
    -o${OUTPUT_PREFIX}/chained/${nb}_vs_${na}.chained.1aln \
    ${OUTPUT_PREFIX}/alns/${nb}_vs_${na}.1aln

ALNplot -p:${nb}_vs_${na} ${OUTPUT_PREFIX}/chained/${nb}_vs_${na}.chained.1aln
mv ${nb}_vs_${na}.pdf ${OUTPUT_PREFIX}/plots/

# Target (Jin668) coordinates: PAF columns 6, 8, 9
ALNtoPAF -T32 ${OUTPUT_PREFIX}/chained/${nb}_vs_${na}.chained.1aln \
    | tee ${OUTPUT_PREFIX}/pafs/${nb}_vs_${na}.paf \
    | awk 'BEGIN{OFS="\t"} {print $6, $8, $9}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > ${OUTPUT_PREFIX}/beds/${nb}_vs_${na}.bed
