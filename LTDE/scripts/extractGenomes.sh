#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=30gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

cd /N/dc2/projects/muri2/Task2/LTDE/data

gunzip -c SoilGenomes_Assembly20151005.tar.gz | tar xvf -
