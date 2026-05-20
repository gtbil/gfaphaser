#!/bin/bash
#SBATCH -N1                     # n processor core(s) per node X 2 threads per core
#SBATCH -n48
#SBATCH --mem=374GB               # maximum memory per node
#SBATCH -p atlas                  # standard node(s)
#SBATCH -J finish_gfa
#SBATCH --time 24:00:00
# INHERIT ASM

source /home/${USER}/.bashrc

cd subgraphs/${ASM}/gfa_walks

vg combine ${ASM}.*.vg > combined.vg
vg paths -E -x  combined.vg > combined.vg.paths

awk '{if ($2 > 100000) {print $1}}' combined.vg.paths > combined.vg.paths.pass
vg paths --retain-paths --paths-file combined.vg.paths.pass -x combined.vg > combined.1.vg
vg mod --remove-non-path combined.1.vg > combined.2.vg

vg convert --gfa-out combined.2.vg > combined.gfa
vg paths --extract-fasta -x combined.2.vg > combined.fa
samtools faidx combined.fa

samtools faidx --region-file <(cut -f1 combined.fa.fai | fgrep "h1_component" | sort) --output hap1.fa combined.fa 
samtools faidx --region-file <(cut -f1 combined.fa.fai | fgrep "h2_component" | sort) --output hap2.fa combined.fa 
