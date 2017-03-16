#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=48:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

#module load samtools
module load python

#mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials
#mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/D100

#for sample in /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/D100/*;
#do
#  line="$(echo "$sample" | cut -d "/" -f11-11)"
#  OUT="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/D100/${line}"
#  #mkdir -p $OUT
#  #cp "${sample}/output/evidence/annotated.gd" "${OUT}/annotated.gd"
#  #cp "${sample}/output/evidence/evidence.gd" "${OUT}/evidence.gd"
#  #cp "${sample}/output/output.gd" "${OUT}/output.gd"
#  #samtools mpileup "${sample}/data/reference.bam" > "${sample}/data/coverage.txt"
#  #python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/getCoverage.py -c -i \
#  #  "${sample}/data/coverage.txt" -o "${sample}/data/coverage_merged.txt"
#  #"${OUT}/evidence.gd"
#done

#echo "Script ran"
python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/calculateNum.py
#python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/test.py
