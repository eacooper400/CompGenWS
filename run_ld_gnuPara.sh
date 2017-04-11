#!/bin/bash
#PBS -N LD
#PBS -l select=1:ncpus=8:mem=15gb,walltime=1:00:00
#PBS -j oe

echo "START ------------------------------"

module load R/3.3.2
module load gnu-parallel

src=/home/ecoope4/CompGenWS

cp $src/Sorghum_c1_r[1-8].phase.vcf /local_scratch/

parallel -j 8 R --vanilla --slave \
	 --args /local_scratch/{} {.}.ld.txt \
	 <$src/12_LD_cluster.R ::: *phase.vcf

cp /local_scratch/*ld.txt $src/

echo "FINISH ----------------------------"

