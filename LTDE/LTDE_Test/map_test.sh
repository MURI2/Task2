#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=10gb,walltime=10:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load cutadapt
module load bwa/0.7.2
module load samtools/0.1.19
module load vcftools


#cd /N/dc2/projects/muri2/Task2/reference_assemblies/ \
#    Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/

bwa index /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta

cd /N/dc2/projects/muri2/Task2/LTDE_Test

# map each set of reads

bwa mem -t 4 /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta \
    ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
    > ./map_output/GSF1018-Lennon_S57_001_NOadapt.sam

bwa mem -t 4 /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta \
    ./GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz ./GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz \
    > ./map_output/GSF1018-Lennon_S97_001_NOadapt.sam

# SAM to BAM
samtools view -Sb  ./map_output/GSF1018-Lennon_S57_001_NOadapt.sam \
    >  ./map_output/GSF1018-Lennon_S57_001_NOadapt.bam

samtools view -Sb ./map_output/GSF1018-Lennon_S97_001_NOadapt.sam \
    >  ./map_output/GSF1018-Lennon_S97_001_NOadapt.bam

cd ./map_output
#samtools merge GSF1018_merged.bam \
#    GSF1018-Lennon_S97_001_NOadapt.bam GSF1018-Lennon_S57_001_NOadapt.bam

samtools merge GSF1018_merged.bam GSF1018-Lennon_S57_001_NOadapt.bam GSF1018-Lennon_S97_001_NOadapt.bam

# Sort the BAM file using SAMtools
samtools sort -o ./GSF1018_merged_sorted.bam ./GSF1018_merged.bam
# Remove PCR duplicates
samtools rmdup ./GSF1018_merged_sorted.bam ./GSF1018_merged_sorted_NOdup.bam
# Index our final BAM file
samtools index ./GSF1018_merged_sorted_NOdup.bam




###### Get the file for mapgd
samtools mpileup -q 25 -Q 25 -B GSF1018_merged_sorted_NOdup.bam \
    > GSF1018_merged_sorted_NOdup.mpileup

samtools view -H GSF1018_merged_sorted_NOdup.bam > GSF1018_merged_sorted_NOdup.header
