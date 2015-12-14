#!/bin/bash
set -e
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
#module load java
#module load fastqc
#module load bioperl
#module load python
#module load gcc
#module load cutadapt
#module load khmer
#module load spades

#PATH=$PATH:/N/soft/mason/galaxy-apps/fastx_toolkit_0.0.13

# set our defaults for the options
KVALUE=51
CUTOFF=auto
# parse the options
while getopts 'c:k:' opt ; do
  case $opt in
    c) CUTOFF=$OPTARG ;;
    k) KVALUE=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1))
# check for mandatory positional parameters
if [ $# -lt 3 ]; then
  echo "Usage: $0 [options] outputFolder Read1 Read2"
  echo "Options: -k kmervalue (def: $KVALUE) | -c cov_cutoff (def: $CUTOFF)"
  exit 1
fi
DIR="$1"
LEFT="$2"
RIGHT="$3"
# do the assembly
echo "Assembling with K=$KVALUE and cutoff=$CUTOFF"
velveth "$DIR" "$KVALUE" -shortPaired -fmtAuto -separate "$LEFT" "$RIGHT"
velvetg "$DIR" -exp_cov auto -cov_cutoff "$CUTOFF" -very_clean yes
echo "Results are in: $DIR"
