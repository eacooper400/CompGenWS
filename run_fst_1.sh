#!/bin/bash
#PBS -N Fst
#PBS -l select=1:ncpus=1:mem=11gb,walltime=0:30:00
#PBS -j oe

echo "START ------------------------------"

module load R/3.3.2
src=/home/ecoope4/CompGenWS

cp $src/4_sampleData_Fst.vcf /local_scratch/

R --vanilla --slave --args /local_scratch/4_sampleData_Fst.vcf /local_scratch/Fst_cluster.txt <$src/12_fst_cluster.R

cp /local_scratch/Fst_cluster.txt $src/

echo "FINISH ----------------------------"

