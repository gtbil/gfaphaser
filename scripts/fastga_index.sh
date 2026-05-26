#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=96G
#SBATCH -t 6:00:00
#SBATCH -J fastga_index
#SBATCH -A gbru_malvaceae
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

# Submit with:
#   GENOMES=mydir sbatch --export=ALL --array=1-$N fastga_index.sh
# Defaults: GENOMES=genomes
#
# Setup before submitting:
#   GENOMES=${GENOMES:-genomes}
#   mkdir -p logs
#   ls $GENOMES/ | grep '\.fa$' > genome_list.txt
#   N=$(wc -l < genome_list.txt)
#   GENOMES=$GENOMES sbatch --export=ALL \
#       -o "logs/fastga_index_%A_%a.out" \
#       -e "logs/fastga_index_%A_%a.err" \
#       --array=1-$N fastga_index.sh
#
# Note: Slurm directives are read at submit time and don't see env vars, so
# -o/-e are passed on the sbatch command line if you want them under logs/.

set -euo pipefail

GENOMES="${GENOMES:-genomes}"

fa="${GENOMES}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" genome_list.txt)"
echo "[task ${SLURM_ARRAY_TASK_ID}] processing ${fa}"

FAtoGDB "$fa"
GIXmake -T12 -P"${TMPDIR}" "${fa%.fa}"
