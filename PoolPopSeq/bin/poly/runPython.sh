#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=8:00:00
#PBS -M wrshoema@umail.iu.edu
#PBS -m abe
#PBS -j oe

module load python

python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/test.py
