#!/bin/bash
#SBATCH --cpus-per-task 3                     # n processor core(s) per node X 2 threads per core
#SBATCH --mem=24GB               # maximum memory per node
#SBATCH -p atlas                  # standard node(s)
#SBATCH -J get_contigs
#SBATCH --time 8:00:00
#SBATCH --array=1-490

# INHERIT ASM
COMPONENT=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")

source /home/${USER}/.bashrc

mkdir -p work/${ASM}/${COMPONENT}
cd work/${ASM}/${COMPONENT}
cp ../../../subgraphs/${ASM}/gfa/${ASM}.component_${COMPONENT}.gfa step01.gfa

/home/${USER}/.pyenv/shims/python ../../../prune_to_spine.py step01.gfa step02.gfa --verbose
[ ! -s step01.gfa ] && exit 1

vg convert --gfa-in step02.gfa --packed-out > step02.vg

vg mod \
	--compact-ids \
	--normalize --until-normal 100 \
	--simplify \
	--orient-forward \
	step02.vg \
	| vg mod --unchop - \
	> step03.vg

vg convert --gfa-out step03.vg > step03.gfa
[ ! -s step03.gfa ] && exit 1


/home/${USER}/.pyenv/shims/python ../../../prune_ports.py step03.gfa step04.gfa --verbose
[ ! -s step04.gfa ] && exit 1

/home/${USER}/.pyenv/shims/python ../../../prune_to_spine.py step04.gfa step05.gfa
[ ! -s step05.gfa ] && exit 1

vg convert --gfa-in step05.gfa --packed-out > step05.vg

vg mod \
	--compact-ids \
	--normalize --until-normal 100 \
	--simplify \
	--orient-forward \
	--prune-subgraphs --length 50000 \
	step05.vg \
	| vg mod --unchop - \
	> step06.vg

vg convert --gfa-out step06.vg > step06.gfa
[ ! -s step06.gfa ] && exit 1

/home/${USER}/.pyenv/shims/python ../../../prune_ports.py step06.gfa step07.gfa --verbose
[ ! -s step07.gfa ] && exit 1

/home/${USER}/.pyenv/shims/python ../../../prune_to_spine.py step07.gfa step08.gfa --verbose
[ ! -s step08.gfa ] && exit 1

/home/${USER}/.pyenv/shims/python ../../../prune_ports.py step08.gfa step09.gfa --verbose
[ ! -s step09.gfa ] && exit 1

../../../gfa_haps step09.gfa
[ ! -s step09.with_paths.gfa ] && exit 1


# repair the path names
sed -i "s/sample/${ASM}/g" step09.with_paths.gfa
sed -i "s/component/component\_${COMPONENT}/g" step09.with_paths.gfa

vg convert --gfa-in step09.with_paths.gfa --packed-out > ../../../subgraphs/${ASM}/gfa_walks/${ASM}.${COMPONENT}.vg
