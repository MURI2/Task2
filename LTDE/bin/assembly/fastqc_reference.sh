#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load fastqc


for file in /N/dc2/projects/muri2/Task2/D400-100/*/*fastq.gz
do
  dType="$(echo "$file" | cut -d "_" -f 2-2 | cut -d "/" -f 2-10)"
  mkdir -p "/N/dc2/projects/muri2/Task2/LTDE/data/reads_quality/D400_100_quality/${dType}"
  fastqc "$file" --outdir="/N/dc2/projects/muri2/Task2/LTDE/data/reads_quality/D400_100_quality/${dType}"
done
