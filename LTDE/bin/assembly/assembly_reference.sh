#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load bioperl
module load python
module load gcc
#module load khmer
module load velvet
module load VelvetOptimiser
module load spades


for folder in /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/*/ ;
do
  R1="${folder}"*"R1_001_cleaned.fastq.gz"
  R2="${folder}"*"R2_001_cleaned.fastq.gz"
  strain="$(echo $folder | cut -d "/" -f 10-10 | cut -d "v" -f 2-2)"
  if [ "$strain" = "13483" ]
  then
    strain='ATCC13985'
  elif [ "$strain" = "43828" ]
  then
    strain='ATCC43928'
  else
    kbs='KBS0'
    strain=$kbs$strain
  fi
  outdir="/N/dc2/projects/muri2/Task2/LTDE/data/2016_KBSGenomes/${strain}"
  mkdir -p $outdir
  #VelvetOptimiser.pl -s 31 -e 71 -f '-shortPaired -fastq.gz' ${R1} '-shortPaired2 -fastq.gz' ${R2} \
  #  -t 8 -k 'n50*ncon' --p $strain'assembled_VO' --d $outdir
  spades.py --careful -o $outdir --pe1-1 $R1 --pe1-2 $R2
done

#VelvetOptimiser.pl -s 31 -e 71 -f '-shortPaired -fastq ./'$GENOME'.trim.interleaved.trim.fastq.pe.keep.pe' -t 8 -k 'n50*ncon' --p $GENOME'assembled'
