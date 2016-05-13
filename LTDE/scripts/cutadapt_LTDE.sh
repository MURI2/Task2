#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load python
module load gcc
module load cutadapt

cd /N/dc2/projects/muri2/Task2/LTDE/data/reads_raw

#declare -a NAMES=("GSF1046-KBS0703-C_S37_R2_001.fastq.gz" "GSF1046-KBS0703-C_S37_R1_001.fastq.gz" "GSF1047-KBS0703-C_S37_R1_001.fastq.gz")

declare -a ARRAY=()

for i in *fastq.gz
#for i in "${NAMES[@]}"
do
  # Get the part of the filename that's in both R1 & R2 reads
  iType="$(echo "$i" | cut -d "_" -f1-2)"
  #ARRAY=(${ARRAY[@]} $iType)
  #ARRAY+=($iType)
  ARRAY=("${ARRAY[@]}" "$iType")
done

# Remove the duplicates
REMDUP=($(printf "%s\n" "${ARRAY[@]}" | sort | uniq -c | sort -rnk1 | awk '{ print $2 }'))

# run cutadapt on each
for j in "${REMDUP[@]}"
do
  InR1="./${j}_R1_001.fastq.gz"
  InR2="./${j}_R2_001.fastq.gz"
  OutR1="../reads_clean/${j}_R1_001_cleaned.fastq.gz"
  OutR2="../reads_clean/${j}_R2_001_cleaned.fastq.gz"
  cutadapt -q 30 -u 10 -u -10 -o $OutR1 -p $OutR2 $InR1 $InR2
done
