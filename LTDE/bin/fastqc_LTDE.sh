#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=4:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load fastqc

cd /N/dc2/projects/muri2/Task2/LTDE/data/reads_raw

for file in ./*fastq.gz
do
  fastqc "$file" --outdir=../reads_quality
done
