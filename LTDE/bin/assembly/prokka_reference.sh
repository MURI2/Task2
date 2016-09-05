####Note#####
# This script runs on my personal mac. Prokka is not installed on dc2


#scp -r wrshoema@mason.indiana.edu:/N/dc2/projects/muri2/Task2/LTDE/data/2016_KBSGenomes/ /Users/WRShoemaker/Desktop/

mkdir -p /Users/WRShoemaker/Desktop/2016_KBSGenomes_Annotate


for i in /Users/WRShoemaker/Desktop/2016_KBSGenomes/*/ ;
do
  strain="$(echo $i | cut -d "/" -f 6-6 )"
  contigs="/Users/WRShoemaker/Desktop/2016_KBSGenomes/${strain}/contigs.fasta"
  OUT="/Users/WRShoemaker/Desktop/2016_KBSGenomes_Annotate/${strain}"
  mkdir -p $OUT
  #iType="$(echo "$i" | cut -d "_" -f1-1)"
  prokka --compliant --centre I --outdir $OUT --locustag G --prefix G-Chr1 $contigs --force
done

#scp -r /Users/WRShoemaker/Desktop/2016_KBSGenomes_Annotate wrshoema@mason.indiana.edu:/N/dc2/projects/muri2/Task2/LTDE/data
