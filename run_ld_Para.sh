#!/bin/bash
#PBS -N LD
#PBS -l select=1:ncpus=8:mem=15gb,walltime=1:00:00
#PBS -j oe

echo "START ------------------------------"

module load R/3.3.2

src=/home/ecoope4/CompGenWS

for file in $src/Sorghum_c1_r[1-8].phase.vcf
do
    export filename=`basename $file .phase.vcf`
    export outname="$filename.ld.txt"
    R --vanilla --slave --args $file $src/$outname <$src/12_LD_cluster.R &
done

wait

echo "FINISH ----------------------------"

