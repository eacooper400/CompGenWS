#!/bin/bash
#PBS -N Fst
#PBS -l select=1:ncpus=1:mem=11gb,walltime=0:30:00
#PBS -j oe

echo "START ------------------------------"

module load R/3.3.2
src=/home/ecoope4/CompGenWS

R --vanilla --slave --args $src/4_sampleData_Fst.vcf $src/Fst_plot.pdf <$src/12_fst_cluster.R

echo "FINISH ----------------------------"

