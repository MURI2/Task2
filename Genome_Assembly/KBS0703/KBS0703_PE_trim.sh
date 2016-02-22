#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load bioperl
module load python
module load gcc
module load cutadapt


cd /N/dc2/projects/muri2/Task2/GSF-911/KBS0703

cutadapt -q 30 -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTCACGCAATCTCGTAT \
    --minimum-length 20 \
    -u 15 \
    -u -20 \
    -o tmp.1.fastq \
    -p tmp.2.fastq \
    GSF911-Ar_S10_L001_R1_001.fastq.gz GSF911-Ar_S10_L001_R2_001.fastq.gz

cutadapt -q 30 -b GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG \
    --minimum-length 20 \
    -u 15 \
    -u -20 \
    -o GSF911-Ar_S10_L001_R2_001_Q30_U15_UN20.fastq.gz \
    -p GSF911-Ar_S10_L001_R1_001_Q30_U15_UN20.fastq.gz \
    tmp.2.fastq tmp.1.fastq

rm tmp.1.fastq tmp.2.fastq
