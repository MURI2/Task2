#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=48:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load java

cd /N/dc2/projects/muri2/Task2/LTDE/data/map_results/

#REF=/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate/KBS0710/G-Chr1

#IN=/N/dc2/projects/muri2/Task2/LTDE/data/map_results/KBS0710/GSF1046-KBS0710-A_S52_mapped_sort_NOdup_sort

##### change output so you have a sub directory with GATK output


for d in */ ;
do
  dType="$(echo "$d" | cut -d "/" -f1-1)"
  REF="/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate/${dType}/G-Chr1"
  for file in $dType/*;
  do
    if [[ $file == *_mapped_sort_NOdup_sort.bam ]]; then
      NoExt="$(echo "${file%.*}")"
      # parse file

      java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
          /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar \
          I="${NoExt}.bam" \
          O="${NoExt}_fixed.bam" \
          SORT_ORDER=coordinate RGID=Rpal RGLB=bar RGPL=illumina RGSM=test_line \
          RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

      java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
          -jar /N/soft/rhel6/picard/picard-tools-1.107/CreateSequenceDictionary.jar \
          R="${REF}.fna" \
          O="${REF}.dict"

      java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
          /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
          -T RealignerTargetCreator \
          -R "${REF}.fna" \
          -I "${NoExt}_fixed.bam" \
          -o "${NoExt}_fixed.intervals"

      java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
          -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
          -T IndelRealigner \
          -R "${REF}.fna" \
          -I "${NoExt}_fixed.bam" \
          -targetIntervals "${NoExt}_fixed.intervals" \
          --filter_bases_not_stored \
          -o "${NoExt}_fixed_realigned.bam"

      java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
          -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
          -T UnifiedGenotyper \
          -R "${REF}.fna" \
          -I "${NoExt}_fixed_realigned.bam" \
          -glm BOTH -rf BadCigar \
          -o "${NoExt}_fixed_realigned.vcf"

      java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
          -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
          -T BaseRecalibrator -I "${NoExt}_fixed_realigned.bam" \
          -R "${REF}.fna" \
          -rf BadCigar --filter_bases_not_stored -knownSites "${NoExt}_fixed_realigned.vcf" \
          -o "${NoExt}_fixed_realigned.recal_data.grp"


      java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
          -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
          -T PrintReads -rf BadCigar \
          -R "${REF}.fna" \
          -I "${NoExt}_fixed_realigned.bam" \
          -o "${NoExt}_fixed_realigned_mapped.bam" \
          -BQSR "${NoExt}_fixed_realigned.recal_data.grp"

      java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
          -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
          -T UnifiedGenotyper \
          -R "${REF}.fna" \
          -I "${NoExt}_fixed_realigned_mapped.bam"\
          -rf BadCigar \
          -o "${NoExt}_fixed_realigned_mapped.vcf"

      subDir="$(echo "$file" | cut -d "_" -f1-1)"
      mkdir -p $subDir
      cp -v $subDir* $subDir/
    else
      continue
    fi
  done
done
