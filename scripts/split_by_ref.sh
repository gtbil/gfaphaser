#!/bin/bash
#SBATCH -N1                     # n processor core(s) per node X 2 threads per core
#SBATCH -n48
#SBATCH --mem=374GB               # maximum memory per node
#SBATCH -p atlas                  # standard node(s)
#SBATCH -J split
#SBATCH --time 24:00:00
# INHERIT ASM

source /home/${USER}/.bashrc

mkdir -p ./subgraphs/${ASM}/miniasm/

mashmap -r ref.fa -q subgraphs/${ASM}/gfa_walks/hap1.fa \
	-s 50000 --pi 90 -t 46 \
	--filter_mode map --dense \
	--hgFilterAniDiff 1.0 \
	-o ./subgraphs/${ASM}/miniasm/hap1.paf

mashmap -r ref.fa -q subgraphs/${ASM}/gfa_walks/hap2.fa \
        -s 50000 --pi 90 -t 46 \
	--filter_mode map --dense \
	--hgFilterAniDiff 1.0 \
        -o ./subgraphs/${ASM}/miniasm/hap2.paf

./assign_contigs.awk ./subgraphs/${ASM}/miniasm/hap1.paf > ./subgraphs/${ASM}/miniasm/hap1.assign.tsv
./assign_contigs.awk ./subgraphs/${ASM}/miniasm/hap2.paf > ./subgraphs/${ASM}/miniasm/hap2.assign.tsv

for CHR in A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13 PT MT; do
	mkdir -p ./subgraphs/${ASM}/miniasm/${CHR};
	awk -v CHR=${CHR} '{if ($3==CHR) {print $1}}' ./subgraphs/${ASM}/miniasm/hap1.assign.tsv | sort > ./subgraphs/${ASM}/miniasm/${CHR}/hap1.paths;
	awk -v CHR=${CHR} '{if ($3==CHR) {print $1}}' ./subgraphs/${ASM}/miniasm/hap2.assign.tsv | sort > ./subgraphs/${ASM}/miniasm/${CHR}/hap2.paths;
	vg paths --retain-paths --paths-file ./subgraphs/${ASM}/miniasm/${CHR}/hap1.paths -x ./subgraphs/${ASM}/gfa_walks/combined.2.vg > ./subgraphs/${ASM}/miniasm/${CHR}/hap1.tmp.vg;
	vg paths --retain-paths --paths-file ./subgraphs/${ASM}/miniasm/${CHR}/hap2.paths -x ./subgraphs/${ASM}/gfa_walks/combined.2.vg > ./subgraphs/${ASM}/miniasm/${CHR}/hap2.tmp.vg;
	vg paths --extract-fasta -x ./subgraphs/${ASM}/miniasm/${CHR}/hap1.tmp.vg  > ./subgraphs/${ASM}/miniasm/${CHR}/hap1.fa;
	vg paths --extract-fasta -x ./subgraphs/${ASM}/miniasm/${CHR}/hap2.tmp.vg  > ./subgraphs/${ASM}/miniasm/${CHR}/hap2.fa;
	rm ./subgraphs/${ASM}/miniasm/${CHR}/hap1.tmp.vg ./subgraphs/${ASM}/miniasm/${CHR}/hap2.tmp.vg;
	minimap2 -x asm5 -X -t 32 ./subgraphs/${ASM}/miniasm/${CHR}/hap1.fa ./subgraphs/${ASM}/miniasm/${CHR}/hap1.fa > ./subgraphs/${ASM}/miniasm/${CHR}/hap1.paf;
	minimap2 -x asm5 -X -t 32 ./subgraphs/${ASM}/miniasm/${CHR}/hap2.fa ./subgraphs/${ASM}/miniasm/${CHR}/hap2.fa > ./subgraphs/${ASM}/miniasm/${CHR}/hap2.paf;
	miniasm -1 -2 -h 100000 -c1 -e1 -f ./subgraphs/${ASM}/miniasm/${CHR}/hap1.fa ./subgraphs/${ASM}/miniasm/${CHR}/hap1.paf > ./subgraphs/${ASM}/miniasm/${CHR}/hap1.gfa;
	miniasm -1 -2 -h 100000 -c1 -e1 -f ./subgraphs/${ASM}/miniasm/${CHR}/hap2.fa ./subgraphs/${ASM}/miniasm/${CHR}/hap2.paf > ./subgraphs/${ASM}/miniasm/${CHR}/hap2.gfa;
done

cd ./subgraphs/${ASM}/miniasm/

for CHR in A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13; do
	cd ${CHR};
	echo ${CHR};
	# post-process to extract all contigs for each chromosome, using a combination of `S` lines from the gfa and unused entries from the fasta
	samtools faidx hap1.fa;
	grep ^S hap1.gfa | awk '{print ">"$2"\n"$3}' | sed "s/^>/>${CHR}.h1./" | fold -w60 > hap1.combined.fa;
	comm -2 -3 <(cut -f1 hap1.fa.fai | sort) <(grep ^a hap1.gfa | cut -f4 | sort) > hap1.regions;
	samtools faidx --region-file hap1.regions --length 60 hap1.fa | sed "s/^>/>${CHR}.h1./" >> hap1.combined.fa;
	rm hap1.regions;
	samtools faidx hap2.fa;
	grep ^S hap2.gfa | awk '{print ">"$2"\n"$3}' | sed "s/^>/>${CHR}.h2./" | fold -w60 > hap2.combined.fa;
	comm -2 -3 <(cut -f1 hap2.fa.fai | sort) <(grep ^a hap2.gfa | cut -f4 | sort) > hap2.regions;
	samtools faidx --region-file hap2.regions --length 60 hap2.fa | sed "s/^>/>${CHR}.h2./" >> hap2.combined.fa;
	rm hap2.regions;
	cd ..;
done;

# combine fastas!

rm -rf hap1.fa hap2.fa
touch hap1.fa
touch hap2.fa

for CHR in A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13; do
	echo ${CHR};
	cat ${CHR}/hap1.combined.fa >> hap1.fa;
	cat ${CHR}/hap2.combined.fa >> hap2.fa;
done;

samtools faidx hap1.fa
samtools faidx hap2.fa
