#!/bin/bash
# submit_pipeline.sh — chain the four pipeline steps with SLURM dependencies
# Usage: ./submit_pipeline.sh <ASM> <PROJECT>

set -euo pipefail

ASM=${1:?Usage: $0 <ASM> <PROJECT>}
PROJECT=${2:?Usage: $0 <ASM> <PROJECT>}

mkdir -p logs/out

# Step 1: prep the gfa
JID1=$(sbatch --parsable \
    --account=${PROJECT} \
    --export=ASM=${ASM} \
    --output=logs/out/${ASM}_prep_gfa_%j.out \
    --error=logs/out/${ASM}_prep_gfa_%j.err \
    prep_gfa.sh)
echo "Submitted prep_gfa.sh           : $JID1"

# Step 2: make haplotypes (array job, waits for step 1)
# NOTE: the array size depends on files created by step 1, so we can't compute
# it now. We submit a tiny helper that runs *after* step 1 and submits step 2
# with the correct array size, then chains steps 3 and 4 off of it.
JID2=$(sbatch --parsable \
    --dependency=afterok:${JID1} \
    --account=${PROJECT} \
    --export=ASM=${ASM},PROJECT=${PROJECT} \
    --job-name=${ASM}_submit_rest \
    --output=logs/out/${ASM}_submit_rest_%j.out \
    --error=logs/out/${ASM}_submit_rest_%j.err \
    --wrap="bash submit_rest.sh")
echo "Submitted dependent submitter   : $JID2 (waits on $JID1)"

echo "Pipeline submitted for ASM=${ASM}"
