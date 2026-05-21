#!/bin/bash
# submit_rest.sh — submits steps 2, 3, 4 after prep_gfa.sh has finished
# ASM and PROJECT come in via --export from submit_pipeline.sh

set -euo pipefail

: "${ASM:?ASM not set}"
: "${PROJECT:?PROJECT not set}"

NCOMP=$(ls subgraphs/${ASM}/gfa/${ASM}.component* | wc -l)
if [[ "$NCOMP" -lt 1 ]]; then
    echo "ERROR: no component files found for ${ASM}" >&2
    exit 1
fi
echo "Found ${NCOMP} component files for ${ASM}"

# Step 2: haplotype paths per subgraph
JID2=$(sbatch --parsable \
    --export=ASM=${ASM} \
    --array=1-${NCOMP} \
    --account=${PROJECT} \
    --output=logs/out/${ASM}_make_haplotypes_%A_%a.out \
    --error=logs/out/${ASM}_make_haplotypes_%A_%a.err \
    make_haplotypes.sh)
echo "Submitted make_haplotypes.sh    : $JID2"

# Step 3: extract haplotype sequences (waits for whole array)
JID3=$(sbatch --parsable \
    --dependency=afterok:${JID2} \
    --export=ASM=${ASM} \
    --account=${PROJECT} \
    --output=logs/out/${ASM}_finish_gfa_%j.out \
    --error=logs/out/${ASM}_finish_gfa_%j.err \
    finish_gfa.sh)
echo "Submitted finish_gfa.sh         : $JID3 (waits on $JID2)"

# Step 4: split components by reference chromosome
JID4=$(sbatch --parsable \
    --dependency=afterok:${JID3} \
    --export=ASM=${ASM} \
    --account=${PROJECT} \
    --output=logs/out/${ASM}_split_by_ref_%j.out \
    --error=logs/out/${ASM}_split_by_ref_%j.err \
    split_by_ref.sh)
echo "Submitted split_by_ref.sh       : $JID4 (waits on $JID3)"
