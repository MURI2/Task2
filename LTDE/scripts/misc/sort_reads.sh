#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

#cd /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean

declare -a ARRAY=()

for i in *fastq.gz
#for i in "${NAMES[@]}"
do
  # Get the part of the filename that's jus the specie
  iType="$(echo "$i" | cut -d "-" -f1-2)"
  #ARRAY=(${ARRAY[@]} $iType)
  #ARRAY+=($iType)
  ARRAY=("${ARRAY[@]}" "$iType")
done

#echo $ARRAY

REMDUP=($(printf "%s\n" "${ARRAY[@]}" | sort | uniq -c | sort -rnk1 | awk '{ print $2 }'))

# make a directory for each species
for j in "${REMDUP[@]}"
do
  echo $j
  mkdir -p $j
done
# move all cleaned files to species directory
for i in *fastq.gz
#for i in "${NAMES[@]}"
do
  iType="$(echo "$i" | cut -d "-" -f1-2)"
  mv $i $iType/$i
  #grep -L -Z -r $i . | xargs -0 -I{} mv {} $i
done

# got through each directory and map reads to reference

# first, create dictionary of references
# loop through directories

# map to reference and output in folder

# merge BAM files.
