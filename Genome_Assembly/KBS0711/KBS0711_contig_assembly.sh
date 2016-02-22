#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load fastqc
module load bioperl
module load python
module load gcc
module load cutadapt
module load khmer
module load spades


cd /N/dc2/projects/muri2/Task2

spades.py --trusted-contigs \
    ./reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE/contigs.fasta \
    --pe1-1 ./GSF-911/KBS0711/GSF911-711_S1_L001_R1_001_Q30_U15_UN20.fastq.gz \
    --pe1-2 ./GSF-911/KBS0711/GSF911-711_S1_L001_R2_001_Q30_U15_UN20.fastq.gz \
    -o ./reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged
