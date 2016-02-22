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

PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13
cd /N/dc2/projects/muri2/Task2

cd /GSF966/KBS0703


cutadapt -q 30 -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATG \
    --minimum-length 20 \
    -o tmp.1.fastq \
    -p tmp.2.fastq \
    GSF966-1-Arthro-6k_S1_R1_001.fastq.gz GSF966-1-Arthro-6k_S1_R2_001.fastq.gz

  cutadapt -q 30 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC \
    --minimum-length 20 \
    -o GSF966-1-Arthro-6k_S1_R2_001_Q30_U10.fastq.gz \
    -p GSF966-1-Arthro-6k_S1_R1_001_Q30_U10.fastq.gz \
    tmp.2.fastq tmp.1.fastq
