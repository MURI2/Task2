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
cutadapt -q 25 -a GTCGTAGAAGATCTCG \
    -o ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
    ./GSF1018-Lennon_S57_R1_001.fastq.gz



cutadapt -q 25 -a GTCGTAGAAGATCTCG -A GTCGTAGAAGATCTCG \
    -o ./GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz -p ./GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz \
    ./GSF1018-Lennon_S97_R1_001.fastq.gz ./GSF1018-Lennon_S97_R2_001.fastq.gz


# Trim bases
#cutadapt -u 25 -o ./white_snps/KBS0711_white_R1_NOadapt_First25.fastq \
#    ./white_snps/KBS0711_white_R1_NOadapt.fastq
#cutadapt -u 25 -o ./white_snps/KBS0711_white_R2_NOadapt_First25.fastq \
#    ./white_snps/KBS0711_white_R2_NOadapt.fastq

# Index reference
bwa index ./Janthino_Genomes/LBCO00000000.1.fasta
cd ./white_snps
# Run BWA
bwa mem -t 4 ../Janthino_Genomes/LBCO00000000.1.fasta ./KBS0711_white_R1_NOadapt_First25.fastq \
    ./KBS0711_white_R2_NOadapt_First25.fastq > ./KBS0711_white_adapt_Q30_First25.sam
# Sort our alignment
samtools view -bS ./KBS0711_white_adapt_Q30_First25.sam | samtools sort - ./KBS0711_white_adapt_Q30_First25_sorted
# Remove PCR duplicates
samtools rmdup ./KBS0711_white_adapt_Q30_First25_sorted.bam ./KBS0711_white_adapt_Q30_First25_sorted_NOdup.bam
# Index our final BAM file
samtools index ./KBS0711_white_adapt_Q30_First25_sorted_NOdup.bam
# Run mpileup & bcftools
# The mean coverage is ~ 27, but -D is for max, gunna go with recommended setting of 100x
samtools mpileup -uf ../Janthino_Genomes/LBCO00000000.1.fasta
    ./KBS0711_white_adapt_Q30_First25_sorted_NOdup.bam | bcftools view -bvcg - > ./white_snps.raw.bcf
bcftools view ./white_snps.raw.bcf | vcfutils.pl varFilter -D150 > ./white_snps_D150.vcf
bcftools view ./white_snps.raw.bcf | vcfutils.pl varFilter -D100 > ./white_snps_D100.vcf
# No filter?
bcftools view ./white_snps.raw.bcf | vcfutils.pl varFilter > ./white_snps_noCutoff.vcf
# Not using varFilter
bcftools view ./white_snps.raw.bcf > ./white_snps_noCutoff_novarFilter.vcf

# test alternative command
bcftools view ./white_snps.raw.bcf > ./test_output_type_noFilter.vcf

vcftools --vcf ./test_output_type_noFilter.vcf --out ./white_snps_maxMeanDP100 --max-meanDP 100 --minQ 10
