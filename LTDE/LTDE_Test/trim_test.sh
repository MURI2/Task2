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

REF=/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate/G-Chr1.fna


cd /N/dc2/projects/muri2/Task2/LTDE/data/LTDE_Test
# Remove adaptors and quality trim with a Phred score of 30
cutadapt -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCTCGTAT \
    -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGCGTAT \
    -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGGGGTT \
    -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCTCGTTT \
    -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGCGGAT \
    -o ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
    ./GSF1018-Lennon_S57_R1_001.fastq.gz


cutadapt -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCTCGTAT \
    -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCTCGTAT \
    -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGCGTAT \
    -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGTAGAATCGCGTAT \
    -o ./GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz -p ./GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz \
    ./GSF1018-Lennon_S97_R1_001.fastq.gz ./GSF1018-Lennon_S97_R2_001.fastq.gz

bwa index $REF

-F 4 -bT $REF

bwa mem -t 4 $REF ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
    > ./GSF1018-Lennon_S57_001_NOadapt.sam

bwa mem -t 4 $REF \
    ./GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz \
    ./GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz \
    > ./GSF1018-Lennon_S97_001_NOadapt.sam

# SAM to BAM mapped

samtools faidx $REF

samtools view -F 4 -bT $REF ./GSF1018-Lennon_S57_001_NOadapt.sam \
    >  ./GSF1018-Lennon_S57_001_NOadapt_mapped.bam

samtools view -F 4 -bT $REF ./GSF1018-Lennon_S97_001_NOadapt.sam \
    >  ./GSF1018-Lennon_S97_001_NOadapt_mapped.bam


samtools sort ./GSF1018-Lennon_S57_001_NOadapt_mapped.bam ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort
samtools index ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort.bam
samtools rmdup ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort.bam ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort_NOdup.bam
samtools index ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort_NOdup.bam
samtools sort ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort_NOdup.bam  ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort_NOdup_sort
samtools index ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort_NOdup_sort.bam

samtools sort ./GSF1018-Lennon_S97_001_NOadapt_mapped.bam ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort
samtools index ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort.bam
samtools rmdup ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort.bam ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort_NOdup.bam
samtools index ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort_NOdup.bam
samtools sort ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort_NOdup.bam  ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort_NOdup_sort
samtools index ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort_NOdup_sort.bam

# merge and move a copy
samtools merge ./GSF1018-Lennon_S97_001_NOadapt_mapped_sort_NOdup_sort.bam \
    ./GSF1018-Lennon_S57_001_NOadapt_mapped_sort_NOdup_sort.bam \
    ./GSF1018-Lennon_S97_S57_001_NOadapt_mapped_sort_NOdup_sort.bam -f

cp ./GSF1018-Lennon_S97_S57_001_NOadapt_mapped_sort_NOdup_sort.bam \
    /N/dc2/projects/muri2/Task2/LTDE/data/map_results/KBS0711/GSF1018-Lennon_S97_S57_001_NOadapt_mapped_sort_NOdup_sort.bam

# SAM to BAM unmapped
samtools view -f 4 -bT $REF ./GSF1018-Lennon_S97_001_NOadapt.sam \
    >  ./GSF1018-Lennon_S97_001_NOadapt_unmapped.bam

samtools view -f 4 -bT $REF ./GSF1018-Lennon_S97_001_NOadapt.sam \
    >  ./GSF1018-Lennon_S97_001_NOadapt_unmapped.bam

samtools sort ./GSF1018-Lennon_S57_001_NOadapt_unmapped.bam ./GSF1018-Lennon_S57_001_NOadapt_unmapped_sort
samtools index ./GSF1018-Lennon_S57_001_NOadapt_unmapped_sort.bam
samtools rmdup ./GSF1018-Lennon_S57_001_NOadapt_unmapped_sort.bam ./GSF1018-Lennon_S57_001_NOadapt_unmapped_sort_NOdup.bam
samtools index ./GSF1018-Lennon_S57_001_NOadapt_unmapped_sort_NOdup.bam
samtools sort ./GSF1018-Lennon_S57_001_NOadapt_unmapped_sort_NOdup.bam  ./GSF1018-Lennon_S57_001_NOadapt_unmapped_sort_NOdup_sort
samtools index ./GSF1018-Lennon_S57_001_NOadapt_unmapped_sort_NOdup_sort.bam

samtools sort ./GSF1018-Lennon_S97_001_NOadapt_unmapped.bam ./GSF1018-Lennon_S97_001_NOadapt_unmapped_sort
samtools index ./GSF1018-Lennon_S97_001_NOadapt_unmapped_sort.bam
samtools rmdup ./GSF1018-Lennon_S97_001_NOadapt_unmapped_sort.bam ./GSF1018-Lennon_S97_001_NOadapt_unmapped_sort_NOdup.bam
samtools index ./GSF1018-Lennon_S97_001_NOadapt_unmapped_sort_NOdup.bam
samtools sort ./GSF1018-Lennon_S97_001_NOadapt_unmapped_sort_NOdup.bam  ./GSF1018-Lennon_S97_001_NOadapt_unmapped_sort_NOdup_sort
samtools index ./GSF1018-Lennon_S97_001_NOadapt_unmapped_sort_NOdup_sort.bam
