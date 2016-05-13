#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load java

cd /N/dc2/projects/muri2/Task2/LTDE/data/map_results/KBS0711/

REF=/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate/G-Chr1

IN=GSF1018-Lennon_S97_S57_001_NOadapt_mapped_sort_NOdup_sort_fixed

# GATK!
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
    /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar \
    I="${IN}.bam" \
    O="${IN}_fixed.bam" \
    SORT_ORDER=coordinate RGID=Rpal RGLB=bar RGPL=illumina RGSM=test_line \
    RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT


#### Create reference dictionary
java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/picard/picard-tools-1.107/CreateSequenceDictionary.jar \
    R="${REF}.fna" \
    O="${REF}.dict"

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
    /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R "${REF}.fna" \
    -I "${IN}_fixed.bam" \
    -o "${IN}_fixed.intervals"

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R "${REF}.fna" \
    -I "${IN}_fixed.bam" \
    -targetIntervals "${IN}_fixed.intervals" \
    --filter_bases_not_stored \
    -o "${IN}_fixed_realigned.bam"

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
    -T UnifiedGenotyper \
    -R "${REF}.fna" \
    -I "${IN}_fixed_realigned.bam" \
    -glm BOTH -rf BadCigar \
    -o "${IN}_fixed_realigned.vcf"

java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -T BaseRecalibrator -I "${IN}_fixed_realigned.bam" \
    -R "${REF}.fna" \
    -rf BadCigar --filter_bases_not_stored -knownSites "${IN}_fixed_realigned.vcf" \
    -o "${IN}_fixed_realigned.recal_data.grp"


java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
    -T PrintReads -rf BadCigar \
    -R "${REF}.fna" \
    -I "${IN}_fixed_realigned.bam" \
    -o "${IN}_fixed_realigned_mapped.bam" \
    -BQSR "${IN}_fixed_realigned.recal_data.grp"
