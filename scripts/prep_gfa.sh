#!/bin/bash
#SBATCH -N1                     # n processor core(s) per node X 2 threads per core
#SBATCH -n48
#SBATCH --mem=374GB               # maximum memory per node
#SBATCH -p atlas                  # standard node(s)
#SBATCH -J prep_gfa
#SBATCH --time 24:00:00
# INHERIT ASM

source /home/${USER}/.bashrc

./get_blunted --provenance ${ASM}.p_utg.txt --threads 46 --verbose --input_gfa ${ASM}.p_utg.gfa > ${ASM}.gfa

vg convert --gfa-in ${ASM}.gfa --packed-out > ${ASM}.vg

vg mod \
	--compact-ids \
	--break-cycles \
	--normalize --until-normal 100 \
	--simplify \
	--unreverse-edges --orient-forward \
	--unchop \
	--prune-subgraphs --length 50000 \
	${ASM}.vg \
	> ${ASM}.norm.vg

vg convert --gfa-out ${ASM}.norm.vg > ${ASM}.norm.gfa
