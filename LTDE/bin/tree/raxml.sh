raxmlHPC-PTHREADS -T 4 -m GTRCAT -p 12345 -x 12345 -# 300 \
  -s /Users/WRShoemaker/github/Task2/LTDE/data/Tree/deathcurves.clustal.afa.fasta \
  -n T14 \
  -w /Users/WRShoemaker/github/Task2/LTDE/data/Tree

raxmlHPC-PTHREADS -T 4 -m GTRCAT -J MRE -# 100 \
  -z /Users/WRShoemaker/github/Task2/LTDE/data/Tree/RAxML_bootstrap.T14 \
  -n T18 -w /Users/WRShoemaker/github/Task2/LTDE/data/Tree
