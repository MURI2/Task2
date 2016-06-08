#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load cutadapt
module load bwa/0.7.2
module load samtools/0.1.19
module load vcftools
module load java



# GATK!
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
    /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar \
    I=test_merge_sorted_NOdup.bam \
    O=test_merge_sorted_NOdup_fixed.bam \
    SORT_ORDER=coordinate RGID=Rpal RGLB=bar RGPL=illumina RGSM=test_line \
    RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT


###### Not working, removing duplicates above via samtools
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
    /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar \
    INPUT=test_merge_sorted_NOdup_fixed.bam \
    OUTPUT=test_merge_sorted_NOdup_fixed_marked.bam \
    M=test_merge_sorted_NOdup_fixed_marked.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
########


#### Create reference dictionary
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/picard/picard-tools-1.107/CreateSequenceDictionary.jar \
    R= /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
    O= /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.dict

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
    /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
    -I test_merge_sorted_NOdup_fixed.bam \
    -o test_merge_sorted_NOdup_fixed.intervals

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
    -I test_merge_sorted_NOdup_fixed.bam \
    -targetIntervals test_merge_sorted_NOdup_fixed.intervals \
    --filter_bases_not_stored \
    -o test_merge_sorted_NOdup_fixed_realigned.bam

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
    -T UnifiedGenotyper \
    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
    -I test_merge_sorted_NOdup_realigned.bam \
    -glm BOTH -rf BadCigar \
    -o test_merge_sorted_NOdup_realigned.vcf

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -T BaseRecalibrator -I test_merge_sorted_NOdup_realigned.bam \
    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
    -rf BadCigar --filter_bases_not_stored -knownSites test_merge_sorted_NOdup_realigned.vcf \
    -o test_merge_sorted_NOdup_realigned.recal_data.grp

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -T PrintReads -rf BadCigar \
    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
    -I test_merge_sorted_NOdup_realigned.bam \
    -o test_merge_sorted_NOdup_realigned_mapped.bam \
    -BQSR test_merge_sorted_NOdup_realigned.recal_data.grp
