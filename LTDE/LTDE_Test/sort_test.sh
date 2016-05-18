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
module load java
module load gatk

cd /N/dc2/projects/muri2/Task2/LTDE_Test
cd ./map_output


samtools faidx /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta
##########Ignore###############
# Remove PCR duplicates
#samtools rmdup ./test_merge_sorted.bam ./test_merge_sorted_NOdups.bam
# Index our final BAM file
#samtools index ./test_merge_sorted_NOdups.bam
##########Ignore###############



# sort bam file
#samtools sort test_merge.bam test_merge_sorted
samtools index test_merge_sorted.bam

###### Get the file for mapgd
#samtools mpileup -q 25 -Q 25 -B test_merge_sorted.bam \
#    > test_merge_sorted.mpileup

#samtools view -H test_merge_sorted.bam > test_merge_sorted.header


# Create a fasta dictionary
#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/picard/picard-tools-1.107/CreateSequenceDictionary.jar \
#    R= /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta \
#    O= /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.dict


# GATK!
#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
#    /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar \
#    I=test_merge_sorted.bam \
#    O=test_merge_sorted_NOdups_fixed.bam \
#    SORT_ORDER=coordinate RGID=Rpal RGLB=bar RGPL=illumina RGSM=test_line \
#    RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
#    /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar \
#    INPUT=test_merge_sorted_NOdups_fixed.bam \
#    OUTPUT=test_merge_sorted_NOdups_fixed_marked.bam \
#    M=test_merge_sorted_NOdups_fixed_marked.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT


java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
    /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta \
    -I test_merge_sorted_NOdups_fixed_marked.bam -T RealignerTargetCreator \
    -o test_merge_sorted_NOdups_fixed_marked.intervals

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta \
#    -I test_merge_sorted_NOdups_fixed_marked.bam -T IndelRealigner \
#    -targetIntervals test_merge_sorted_NOdups_fixed_marked.intervals --filter_bases_not_stored \
#    -o test_merge_sorted_NOdups_fixed_marked_realigned.bam

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
#    -T UnifiedGenotyper \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta \
#    -I test_merge_sorted_NOdups_fixed_marked_realigned.bam \
#    -glm BOTH -rf BadCigar \
#    -o test_merge_sorted_NOdups_fixed_marked_realigned.vcf

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
#    -T BaseRecalibrator -I test_merge_sorted_NOdups_fixed_marked_realigned.bam \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta \
#    -rf BadCigar --filter_bases_not_stored -knownSites test_merge_sorted_NOdups_fixed_marked_realigned.vcf \
#    -o test_merge_sorted_NOdups_fixed_marked_realigned.recal_data.grp

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged/contigs.fasta \
#    -I test_merge_sorted_NOdups_fixed_marked_realigned.bam -T PrintReads -rf BadCigar \
#    -o test_merge_sorted_NOdups_fixed_marked_realigned_mapped.bam \
#     -BQSR test_merge_sorted_NOdups_fixed_marked_realigned.recal_data.grp

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fna -I Sample_${i}.sorted.fixed.marked.realigned.bam -rf BadCigar -o Sample_${i}.sorted.fixed.marked.realigned.vcf" >> Task3${i}_EcoliAssemble.sh
