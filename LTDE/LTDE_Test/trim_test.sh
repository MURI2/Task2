#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=2gb,walltime=10:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load cutadapt
module load bwa/0.7.2
module load samtools/0.1.19
module load vcftools

# So I'm running the whole thing fresh to make sure it works right.

cd /N/dc2/projects/muri2/Task2/LTDE_Test
# Remove adaptors and quality trim with a Phred score of 30
cutadapt -q 30 -a GTCGTAGAAGATCTCG \
    -o ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
    ./GSF1018-Lennon_S57_R1_001.fastq.gz


#cutadapt -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCTCGTAT \
#    -o ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
#    ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz

#cutadapt -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGCGTAT \
#    -o ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
#    ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz


#cutadapt -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGGGGTT \
#    -o ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
#    ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz

#cutadapt -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCTCGTTT \
#    -o ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
#    ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz

#cutadapt -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGCGGAT \
#    -o ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
#    ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz



cutadapt -q 30 -a GTCGTAGAAGATCTCG -A GTCGTAGAAGATCTCG \
    -o ./GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz -p ./GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz \
    ./GSF1018-Lennon_S97_R1_001.fastq.gz ./GSF1018-Lennon_S97_R2_001.fastq.gz

#cutadapt -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCTCGTAT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGCGTAT \
#    -o ./GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz -p ./GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz \
#    ./GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz ./GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz
