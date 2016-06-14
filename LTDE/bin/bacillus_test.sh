#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load bwa/0.7.2
module load samtools/0.1.19
module load vcftools
module load java
module rm gcc
module load gcc/4.9.2
module load gsl/1.15


cd /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean

declare -a SoilGen=("KBS0812")

declare -a ARRAY=()

REF="/N/dc2/projects/muri2/Task2/LTDE/data/Bacillus_test/AL009126.fa"
#bwa index $REF
#samtools faidx $REF

# pass string and array and check if the array contains the string
containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 1; done
  return 0
}


#for i in *fastq.gz
#for i in "${NAMES[@]}"
#do
  # Get the part of the filename that's in both R1 & R2 reads
#  iType="$(echo "$i" | cut -d "_" -f1-2)"

#  ARRAY=("${ARRAY[@]}" "$iType")
#done

# Remove the duplicates
REMDUP=($(printf "%s\n" "${ARRAY[@]}" | sort | uniq -c | sort -rnk1 | awk '{ print $2 }'))



# run bwa on each bacillus sample
#for j in "${REMDUP[@]}"
#do
#  SPECIES="$(echo "$j" |cut -d"-" -f2 | cut -d"-" -f1)"
#  containsElement "$SPECIES" "${SoilGen[@]}"
#  Result="$(echo $?)"
  # 1 == True
#  if [ "$Result" -eq 1 ]; then
#    echo $j
    # get and index the reference
#    mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/map_results/$SPECIES
#    REF="/N/dc2/projects/muri2/Task2/LTDE/data/Bacillus_test/AL009126.fa"
#    R1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${j}_R1_001_cleaned.fastq.gz"
#    R2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${j}_R2_001_cleaned.fastq.gz"
#    OUT1="/N/dc2/projects/muri2/Task2/LTDE/data/map_results/${SPECIES}"
#    bwa mem -t 4 $REF $R1 $R2 > "${OUT1}/${j}.sam"
    # mapped reads
#    samtools view -F 4 -bT $REF "${OUT1}/${j}.sam" >  "${OUT1}/${j}_mapped.bam"
    # unmapped reads
#    samtools view -f 4 -bT $REF "${OUT1}/${j}.sam" >  "${OUT1}/${j}_unmapped.bam"

#    samtools sort "${OUT1}/${j}_mapped.bam" "${OUT1}/${j}_mapped_sort"
#    samtools index "${OUT1}/${j}_mapped_sort.bam"
#    samtools rmdup "${OUT1}/${j}_mapped_sort.bam" "${OUT1}/${j}_mapped_sort_NOdup.bam"
#    samtools index "${OUT1}/${j}_mapped_sort_NOdup.bam"
#    samtools sort "${OUT1}/${j}_mapped_sort_NOdup.bam" "${OUT1}/${j}_mapped_sort_NOdup_sort"
#    samtools index "${OUT1}/${j}_mapped_sort_NOdup_sort.bam"
    # same thing for unmapped reads
#    samtools sort "${OUT1}/${j}_unmapped.bam" "${OUT1}/${j}_unmapped_sort"
#    samtools index "${OUT1}/${j}_unmapped_sort.bam"
#    samtools rmdup "${OUT1}/${j}_unmapped_sort.bam" "${OUT1}/${j}_unmapped_sort_NOdup.bam"
#    samtools index "${OUT1}/${j}_unmapped_sort_NOdup.bam"
#    samtools sort "${OUT1}/${j}_unmapped_sort_NOdup.bam" "${OUT1}/${j}_unmapped_sort_NOdup_sort"
#    samtools index "${OUT1}/${j}_unmapped_sort_NOdup_sort.bam"
#  else
#    continue
#  fi
#done

cd /N/dc2/projects/muri2/Task2/LTDE/data/map_results/

dType="KBS0812"
#mkdir -p "${dType}_MAPGD"
#samtools merge "${dType}_MAPGD/${dType}_merged.bam" *"_mapped_sort_NOdup_sort.bam"
#samtools view -H "${dType}_MAPGD/${dType}_merged.bam" > "${dType}_MAPGD/${dType}_merged.header"
#samtools mpileup -q 25 -Q 25 -B *"_mapped_sort_NOdup_sort.bam" \
#| /N/dc2/projects/muri2/Task2/LTDE/MAPGD-master/bin/mapgd proview -H "${dType}_MAPGD/${dType}_merged.header" \
#| /N/dc2/projects/muri2/Task2/LTDE/MAPGD-master/bin/mapgd pool -a 22 -o "${dType}_MAPGD/${dType}_merged"


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
        R="${REF}" \
        O="/N/dc2/projects/muri2/Task2/LTDE/data/Bacillus_test/AL009126.dict"

    java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar \
        /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R "${REF}" \
        -I "${NoExt}_fixed.bam" \
        -o "${NoExt}_fixed.intervals"

    java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
        -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R "${REF}" \
        -I "${NoExt}_fixed.bam" \
        -targetIntervals "${NoExt}_fixed.intervals" \
        --filter_bases_not_stored \
        -o "${NoExt}_fixed_realigned.bam"

    java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
        -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
        -T UnifiedGenotyper \
        -R "${REF}" \
        -I "${NoExt}_fixed_realigned.bam" \
        -glm BOTH -rf BadCigar \
        -o "${NoExt}_fixed_realigned.vcf"

    java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
        -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
        -T BaseRecalibrator -I "${NoExt}_fixed_realigned.bam" \
        -R "${REF}" \
        -rf BadCigar --filter_bases_not_stored -knownSites "${NoExt}_fixed_realigned.vcf" \
        -o "${NoExt}_fixed_realigned.recal_data.grp"


    java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
        -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
        -T PrintReads -rf BadCigar \
        -R "${REF}" \
        -I "${NoExt}_fixed_realigned.bam" \
        -o "${NoExt}_fixed_realigned_mapped.bam" \
        -BQSR "${NoExt}_fixed_realigned.recal_data.grp"

    java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
        -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 \
        -T UnifiedGenotyper \
        -R "${REF}" \
        -I "${NoExt}_fixed_realigned_mapped.bam"\
        -rf BadCigar \
        -o "${NoExt}_fixed_realigned_mapped.vcf"

    java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ \
      -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar \
      -R "${REF}" \
      -T VariantsToTable \
      -V "${NoExt}_fixed_realigned_mapped.vcf" \
      -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F AC \
      -o "${NoExt}_fixed_realigned_mapped.txt"

    subDir="$(echo "$file" | cut -d "_" -f1-1)"
    mkdir -p $subDir
    mv -v $subDir* $subDir/
  else
    continue
  fi
done
