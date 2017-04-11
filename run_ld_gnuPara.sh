#!/bin/bash
#PBS -N LD
#PBS -l select=1:ncpus=8:mem=15gb,walltime=1:00:00
#PBS -j oe

echo "START ------------------------------"

module load R/3.3.2
module load gnu-parallel

src=/home/ecoope4/CompGenWS

parallel -j 8 R $src/12_LD_cluster.R {} {.}.ld.txt ::: $src/Sorghum_c1_r1.phase.vcf $src/Sorghum_c1_r2.phase.vcf $src/Sorghum_c1_r3.phase.vcf $src/Sorghum_c1_r4.phase.vcf $src/Sorghum_c1_r5.phase.vcf



echo "FINISH ----------------------------"

