#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load samtools/0.1.19

cd /N/dc2/projects/muri2/Task2/LTDE/data/map_results

for i in */*/*_mapped_sort_NOdup_sort.bam
do
  NoExt="$(echo "${i%.*}")"
  samtools view -h -o "${NoExt}.sam" $i
done
