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

cd /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean

declare -a SoilGen=("KBS0701" "KBS0705" "KBS0710" "KBS0713" "KBS0715" "KBS0722"
"KBS0725" "KBS0802" "KBS0703" "KBS0706" "KBS0711" "KBS0714" "KBS0721" "KBS0724" "KBS0727")

declare -a ARRAY=()

# index the reference genomes that have beeen annotated
for i in "${SoilGen[@]}"
do
  REF="/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate/${i}/G-Chr1.fna"
  bwa index $REF
  samtools faidx $REF
done

# pass string and array and check if the array contains the string
containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 1; done
  return 0
}


for i in *fastq.gz
#for i in "${NAMES[@]}"
do
  # Get the part of the filename that's in both R1 & R2 reads
  iType="$(echo "$i" | cut -d "_" -f1-2)"
  #ARRAY=(${ARRAY[@]} $iType)
  #ARRAY+=($iType)
  ARRAY=("${ARRAY[@]}" "$iType")
done

# Remove the duplicates
REMDUP=($(printf "%s\n" "${ARRAY[@]}" | sort | uniq -c | sort -rnk1 | awk '{ print $2 }'))


# run bwa on each
for j in "${REMDUP[@]}"
do
  SPECIES="$(echo "$j" |cut -d"-" -f2 | cut -d"-" -f1)"
  containsElement "$SPECIES" "${SoilGen[@]}"
  Result="$(echo $?)"
  # 1 == True
  if [ "$Result" -eq 1 ]; then
    echo $j
    # get and index the reference
    mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/map_results/$SPECIES
    REF="/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate/${SPECIES}/G-Chr1.fna"
    R1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${j}_R1_001_cleaned.fastq.gz"
    R2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${j}_R2_001_cleaned.fastq.gz"
    OUT1="/N/dc2/projects/muri2/Task2/LTDE/data/map_results/${SPECIES}"
    bwa mem -t 4 $REF $R1 $R2 > "${OUT1}/${j}.sam"
    # mapped reads
    samtools view -F 4 -bT $REF "${OUT1}/${j}.sam" >  "${OUT1}/${j}_mapped.bam"
    # unmapped reads
    samtools view -f 4 -bT $REF "${OUT1}/${j}.sam" >  "${OUT1}/${j}_unmapped.bam"

    samtools sort "${OUT1}/${j}_mapped.bam" "${OUT1}/${j}_mapped_sort"
    samtools index "${OUT1}/${j}_mapped_sort.bam"
    samtools rmdup "${OUT1}/${j}_mapped_sort.bam" "${OUT1}/${j}_mapped_sort_NOdup.bam"
    samtools index "${OUT1}/${j}_mapped_sort_NOdup.bam"
    samtools sort "${OUT1}/${j}_mapped_sort_NOdup.bam" "${OUT1}/${j}_mapped_sort_NOdup_sort"
    samtools index "${OUT1}/${j}_mapped_sort_NOdup_sort.bam"
    samtools view -h -o "${OUT1}/${j}_mapped_sort_NOdup_sort.sam" "${OUT1}/${j}_mapped_sort_NOdup_sort.bam"
    # same thing for unmapped reads
    samtools sort "${OUT1}/${j}_unmapped.bam" "${OUT1}/${j}_unmapped_sort"
    samtools index "${OUT1}/${j}_unmapped_sort.bam"
    samtools rmdup "${OUT1}/${j}_unmapped_sort.bam" "${OUT1}/${j}_unmapped_sort_NOdup.bam"
    samtools index "${OUT1}/${j}_unmapped_sort_NOdup.bam"
    samtools sort "${OUT1}/${j}_unmapped_sort_NOdup.bam" "${OUT1}/${j}_unmapped_sort_NOdup_sort"
    samtools index "${OUT1}/${j}_unmapped_sort_NOdup_sort.bam"

  else
    continue
  fi
done
