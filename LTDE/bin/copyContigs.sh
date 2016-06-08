#!/bin/bash

declare -a SoilGen=("KBS0701" "KBS0705" "KBS0710" "KBS0713" "KBS0715" "KBS0722"
"KBS0725" "KBS0802" "KBS0703" "KBS0706" "KBS0711" "KBS0714" "KBS0721" "KBS0724" "KBS0727")

declare -a ARRAY=()

# make directory for prokka
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate

# index the reference genomes
for i in "${SoilGen[@]}"
do
  REF="/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes/${i}/${i}assembled_data_"*/"contigs.fa"
  cp $REF "/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Annotate/${i}_contigs.fa"
done
