#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=10gb,walltime=10:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load java
module load fastqc

cd /N/dc2/projects/muri2/Task2/LTDE_Test

fastqc GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz
fastqc GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz
fastqc GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz
