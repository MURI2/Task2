#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load python
module load gcc
module load cutadapt

#-g CAAGCAGAAGACGGCATACGA
#-g AATGATACGGCGACCACCGA

mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean"
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/D100"


for folder in /N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/D100/*/
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
  mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/D100/${line_folder}"

  for rep in "${REMDUPrep[@]}"
  do
    InR1="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/D100/${REMDUPline}_R1_${rep}.fastq.gz"
    InR2="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_raw/D100/${REMDUPline}_R2_${rep}.fastq.gz"
    OutR1="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/D100/${REMDUPline}_R1_${rep}_clean.fastq.gz"
    OutR2="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean/D100/${REMDUPline}_R2_${rep}_clean.fastq.gz"
    adaptor1="$(  echo "$REMDUPline" | cut -d"_" -f3-3 | cut -d"-" -f1-1)"
    adaptor2="$(  echo "$REMDUPline" | cut -d"_" -f3-3 | cut -d"-" -f2-2)"
    #echo $adaptor1
    #echo $adaptor2
    #   -b AGATCGGAAGAGC -B AGATCGGAAGAGC -b $adaptor1 -B $adaptor2 \
    # -b AATGATACGGCGACCACCGA -B CAAGCAGAAGACGGCATACGA \
    #-b "CAAGCAGAAGACGGCATACGAGAT${adaptor1}GTCTCGTGGGCTCGG" \
    #-B "AATGATACGGCGACCACCGAGATCTACAC${adaptor2}TCGTCGGCAGCGTC" \
    #-b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -B GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
    cutadapt -q 30,30 -u 10 \
      -b $adaptor1 -B $adaptor1 -b $adaptor2 -B $adaptor2 \
      -b TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -B GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
      -o $OutR1 -p $OutR2 $InR1 $InR2
  done
done
