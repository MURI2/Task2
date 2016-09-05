####Note#####
# This script runs on my personal mac. Prokka is not installed on dc2

cd /Users/WRShoemaker/Desktop/

scp -r wrshoema@mason.indiana.edu:/N/dc2/projects/muri2/Task2/LTDE/data/2015_SoilGenomes_Assembled/ .

mkdir -p 2015_SoilGenomes_Annotate/2015_SoilGenomes_Annotate

cd 2015_SoilGenomes_Annotate/2015_SoilGenomes_Assembled

for i in *.fa
do
  iType="$(echo "$i" | cut -d "_" -f1-1)"
  OUT="/Users/WRShoemaker/Desktop/2015_SoilGenomes_Annotate/${iType}"
  prokka --compliant --centre I --outdir $OUT --locustag G --prefix G-Chr1 $i --force
done

scp -r /Users/WRShoemaker/Desktop/2015_SoilGenomes_Annotate wrshoema@mason.indiana.edu:/N/dc2/projects/muri2/Task2/LTDE/data
