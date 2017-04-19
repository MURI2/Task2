#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load java
module load python
module load gcc
module load cutadapt

#-g CAAGCAGAAGACGGCATACGA
#-g AATGATACGGCGACCACCGA

Sample_date=D100

mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic"
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${Sample_date}"


for folder in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/${Sample_date}/"*/
do
  declare -a ARRAYreps=()
  declare -a ARRAYlines=()

  for file in $folder*.fastq.gz
  do
    rep="$(  echo "$file" | cut -d'.' -f 1 | rev | cut -d"_" -f1 | rev)"
    ARRAYreps=("${ARRAYreps[@]}" "$rep")
    line="$(  echo "$file" | cut -d"_" -f2-5 | cut -d"/" -f3-5)"
    ARRAYlines=("${ARRAYlines[@]}" "$line")
  done
  REMDUPrep=($(printf "%s\n" "${ARRAYreps[@]}" | sort | uniq -c | sort -rnk1 | awk '{ print $2 }'))
  REMDUPline=($(printf "%s\n" "${ARRAYlines[@]}" | sort | uniq -c | sort -rnk1 | awk '{ print $2 }'))
  line_folder="$(  echo "$REMDUPline" | cut -d"/" -f1-1)"
  mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${Sample_date}/${line_folder}"

  for rep in "${REMDUPrep[@]}"
  do
    InR1="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/${Sample_date}/${REMDUPline}_R1_${rep}.fastq.gz"
    InR2="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/${Sample_date}/${REMDUPline}_R2_${rep}.fastq.gz"
    OutR1Paired="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${Sample_date}/${REMDUPline}_R1_${rep}_clean_paired.fastq.gz"
    OutR2Paired="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${Sample_date}/${REMDUPline}_R2_${rep}_clean_paired.fastq.gz"
    OutR1UnPaired="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${Sample_date}/${REMDUPline}_R1_${rep}_clean_unpaired.fastq.gz"
    OutR2UnPaired="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${Sample_date}/${REMDUPline}_R2_${rep}_clean_unpaired.fastq.gz"
    adaptor1="$(  echo "$REMDUPline" | cut -d"_" -f3-3 | cut -d"-" -f1-1)"
    adaptor2="$(  echo "$REMDUPline" | cut -d"_" -f3-3 | cut -d"-" -f2-2)"
    #echo $adaptor1
    #echo $adaptor2
    java -jar /N/dc2/projects/MicroEukMA/softwares/Trimmomatic-0.32/trimmomatic-0.32.jar \
      PE -threads 4 $InR1 $InR2 $OutR1Paired $OutR1UnPaired $OutR2Paired $OutR2UnPaired \
      ILLUMINACLIP:/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/transposase.fa:2:30:10 \
      LEADING:4 TRAILING:4 MINLEN:40 HEADCROP:15
  done
done
