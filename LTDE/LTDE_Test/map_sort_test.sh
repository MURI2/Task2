#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
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

REF=/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate/G-Chr1.fna

bwa index $REF

cd /N/dc2/projects/muri2/Task2/LTDE_Test

bwa mem -t 4 /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
    ./GSF1018-Lennon_S57_R1_001_NOadapt.fastq.gz \
    > ./map_output/GSF1018-Lennon_S57_001_NOadapt.sam

bwa mem -t 4 /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
    ./GSF1018-Lennon_S97_R1_001_NOadapt.fastq.gz ./GSF1018-Lennon_S97_R2_001_NOadapt.fastq.gz \
    > ./map_output/GSF1018-Lennon_S97_001_NOadapt.sam

# SAM to BAM
samtools view -Sb  ./map_output/GSF1018-Lennon_S57_001_NOadapt.sam \
    >  ./map_output/GSF1018-Lennon_S57_001_NOadapt.bam

samtools view -Sb ./map_output/GSF1018-Lennon_S97_001_NOadapt.sam \
    >  ./map_output/GSF1018-Lennon_S97_001_NOadapt.bam

cd ./map_output

#### THis is the good command
#### -f = force to overwrite the file
samtools merge GSF1018_merged.bam GSF1018-Lennon_S57_001_NOadapt.bam GSF1018-Lennon_S97_001_NOadapt.bam -f
######

samtools sort ./GSF1018_merged.bam ./GSF1018_merged_sorted

samtools index ./GSF1018_merged_sorted.bam


cd /N/dc2/projects/muri2/Task2/LTDE_Test/map_output


samtools faidx /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna


samtools rmdup ./GSF1018_merged_sorted.bam ./GSF1018_merged_sorted_NOdup.bam

#samtools index GSF1018_merged_sorted_NOdup.bam

samtools sort GSF1018_merged_sorted_NOdup.bam ./GSF1018_merged_sorted_NOdup_sorted

samtools index GSF1018_merged_sorted_NOdup_sorted.bam



# GATK!
#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
#    /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar \
#    I=test_merge_sorted_NOdup.bam \
#    O=test_merge_sorted_NOdup_fixed.bam \
#    SORT_ORDER=coordinate RGID=Rpal RGLB=bar RGPL=illumina RGSM=test_line \
#    RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT


###### Not working, removing duplicates above via samtools
#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
#    /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar \
#    INPUT=test_merge_sorted_NOdup_fixed.bam \
#    OUTPUT=test_merge_sorted_NOdup_fixed_marked.bam \
#    M=test_merge_sorted_NOdup_fixed_marked.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
########


#### Create reference dictionary
#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/picard/picard-tools-1.107/CreateSequenceDictionary.jar \
#    R= /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
#    O= /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.dict

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
#    /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
#    -T RealignerTargetCreator \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
#    -I test_merge_sorted_NOdup_fixed.bam \
#    -o test_merge_sorted_NOdup_fixed.intervals

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
#    -T IndelRealigner \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
#    -I test_merge_sorted_NOdup_fixed.bam \
#    -targetIntervals test_merge_sorted_NOdup_fixed.intervals \
#    --filter_bases_not_stored \
#    -o test_merge_sorted_NOdup_fixed_realigned.bam

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
#    -T UnifiedGenotyper \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
#    -I test_merge_sorted_NOdup_realigned.bam \
#    -glm BOTH -rf BadCigar \
#    -o test_merge_sorted_NOdup_realigned.vcf

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
#    -T BaseRecalibrator -I test_merge_sorted_NOdup_realigned.bam \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
#    -rf BadCigar --filter_bases_not_stored -knownSites test_merge_sorted_NOdup_realigned.vcf \
#    -o test_merge_sorted_NOdup_realigned.recal_data.grp

#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
#    -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
#    -T PrintReads -rf BadCigar \
#    -R /N/dc2/projects/muri2/Task2/reference_assemblies/Janthinobacterium_sp_KBS0711/KBS0711_6kMP_11kMP_PE_merged_PROKKA/X_03292016.fna \
#    -I test_merge_sorted_NOdup_realigned.bam \
#    -o test_merge_sorted_NOdup_realigned_mapped.bam \
#    -BQSR test_merge_sorted_NOdup_realigned.recal_data.grp


##### MAPGD ########

module rm gcc
module load gcc/4.9.2
module load gsl/1.15



# try the long way
#samtools mpileup -q 25 -Q 25 -B \
#    /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/test_merge_sorted_NOdup.bam \
#    > /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population.mpileup

#samtools view -H \
#    /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/test_merge_sorted_NOdup.bam \
#    > /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population.header

#bin/mapgd proview -i \
#    /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population.mpileup \
#    -H /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population.header \
#    > /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population.pro

#bin/mapgd pool -i \
#    /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population.pro \
#    -o /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population_allelfrequencies.map



#samtools mpileup -q 25 -Q 25 -B /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/test_merge_sorted_NOdup.bam \
#    | bin/mapgd proview -H /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population.header \
#    | bin/mapgd pool -a 22 -o allelefrequency


##########
#samtools view -H \
#    /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/test_merge_sorted_NOdup.bam \
#    > /N/dc2/projects/muri2/Task2/LTDE_Test/map_output/population.header


#### good !!!!
cd /N/dc2/projects/muri2/Task2/LTDE_Test/map_output

samtools view -H GSF1018_merged_sorted_NOdup_sorted.bam > GSF1018_merged_sorted_NOdup_sorted.header

/N/dc2/projects/muri2/Tools/MAPGD/bin/mapgd sam2idx -H GSF1018_merged_sorted_NOdup_sorted.header > GSF1018_merged_sorted_NOdup_sorted.idx

samtools mpileup -q 25 -Q 25 -B GSF1018_merged_sorted_NOdup_sorted.bam > GSF1018_merged_sorted_NOdup_sorted.mpileup

#/N/dc2/projects/muri2/Tools/MAPGD/bin/mapgd proview \
#    -i test_merge_sorted_NOdup_sorted.mpileup \
#    -H test_merge_sorted_NOdup_sorted.header \
#    > test_merge_sorted_NOdup_sorted.pro

#/N/dc2/projects/muri2/Tools/MAPGD/bin/mapgd pool \
#    -i test_merge_sorted_NOdup_sorted.pro \
#    -o test_merge_sorted_NOdup_sorted.map


#samtools mpileup -q 25 -Q 25 -B test_merge_sorted_NOdup_sorted.bam \
#    | /N/dc2/projects/muri2/Tools/MAPGD/bin/mapgd proview -H test_merge_sorted_NOdup_sorted.idx \
#    | /N/dc2/projects/muri2/Tools/MAPGD/bin/mapgd pool -a 22 -o test_merge_sorted_NOdup_sorted.pol

############ OMGGGGGGG!!!!
/N/dc2/projects/muri2/Tools/MAPGD/bin/mapgd proview \
    -i GSF1018_merged_sorted_NOdup_sorted.mpileup \
    -H GSF1018_merged_sorted_NOdup_sorted.header \
    | /N/dc2/projects/muri2/Tools/MAPGD/bin/mapgd pool -a 22 -o GSF1018_merged_sorted_NOdup_sorted.pol
##### good!!!
