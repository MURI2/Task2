#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=5:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load python
module load gcc
module load cutadapt

for folder in /N/dc2/projects/muri2/Task2/D400-100/*/ ;
do
  fileName="$(echo "${folder}" | cut -d "/" -f 8-8)"
  mkdir -p "/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${fileName}"
  for file in "${folder}"*"R1_001.fastq.gz" ;
  do
    dType="$(echo "${file}" | cut -d "_" -f 1-4)"
    dType1="$(echo "$file" | cut -d "_" -f 2-4 | cut -d "/" -f 2-10)"
    INR1="${dType}_R1_001.fastq.gz"
    INR2="${dType}_R2_001.fastq.gz"
    OutR1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${fileName}/${dType1}_R1_001_cleaned.fastq.gz"
    OutR2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${fileName}/${dType1}_R2_001_cleaned.fastq.gz"
    cutadapt -q 30 -u 15 -o $OutR1 -p $OutR2 $INR1 $INR2
  done
done
