#!/bin/bash

#cd /N/dc2/projects/muri2/Task2


if [ $# -lt 9 ]; then
  echo "Usage: $0 outputFolder PEL PER MP1L MP1R MP2L MP2R PhredScore Cutoff"
  exit 1
fi

OUT="$1"
PEL="$2"
PER="$3"
MP1L="$4"
MP1R="$5"
MP2L="$6"
MP2R="$7"
PHRED="8"
CUTOFF="$9"



#cutadapt -q "$PHRED" -b

echo \
CUTOFF
