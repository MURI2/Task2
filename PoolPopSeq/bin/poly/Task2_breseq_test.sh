#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load bowtie2
module load intel
module load curl
module load java
module load R
module load breseq

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output/D100

OUT=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output/D100/Sample_L0B1
mkdir -p $OUT
REF=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Bacillus_subtilis_168/GCA_000009045.1_ASM904v1_genomic.gbff
SamplePath=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/D100/Sample_L0B1
breseq -j 8 -p -o $OUT -r $REF "${SamplePath}/"*_clean_paired.fastq.gz
