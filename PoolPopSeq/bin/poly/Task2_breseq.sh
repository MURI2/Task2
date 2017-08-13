#!/bin/bash

DAY=D100

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/breseq/line_scripts_gbk
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/${DAY}"
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/breseq/line_scripts_gbk
mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/breseq/line_scripts_gbk/${DAY}"

P=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pseudomonas_sp_KBS0710/G-Chr1.gbk
D=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCA_000008565.1_ASM856v1_genomic.gbff
#B=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Bacillus_subtilis_168/GCA_000009045.1_ASM904v1_genomic.gbff
B=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff
C=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCA_000022005.1_ASM2200v1_genomic.gbff
F=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.gbk
J=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/KBS0711_2015_SoilGenomes_Annotate/G-Chr1.gbk

for line_path in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${DAY}/"*;
do
  line="$(echo "$line_path" | cut -d "/" -f11-11)"
  bash_out="/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/breseq/line_scripts_gbk/${DAY}/${line}_breseq.sh"
  #echo $line
  if [ ! -f $bash_out ]; then
    taxon="$(echo "$line_path" | grep -Po ".(?=.{5}$)")"
    if [[ $taxon == "P" ]]; then
      REF=$P
      #continue
    elif [[ $taxon == "D" ]]
    then
      REF=$D
      #continue
    elif [[ $taxon == "B" ]]
    then
      REF=$B
      #continue
    elif [[ $taxon == "S" ]]
    then
      REF=$B
      #continue
    elif [[ $taxon == "C" ]]
    then
      REF=$C
      #continue
    elif [[ $taxon == "F" ]]
    then
      REF=$F
      #continue
    elif [[ $taxon == "J" ]]
    then
      REF=$J
      #continue
    else
      continue
    fi

    reads="${line_path}/"*_clean_paired.fastq.gz
    OUT="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk/${DAY}/${line}"
    OUT_essentials="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_essentials/${DAY}/${line}"
    mkdir -p $OUT
    mkdir -p $OUT_essentials

    echo '#!/bin/bash' >> $bash_out
    echo '#PBS -k o' >> $bash_out
    echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00' >> $bash_out
    echo '#PBS -M wrshoema@umail.iu.edu' >> $bash_out
    echo '#PBS -m abe' >> $bash_out
    echo '#PBS -j oe' >> $bash_out
    echo '' >> $bash_out
    echo 'module load python' >> $bash_out
    echo 'module load bowtie2/2.2.6' >> $bash_out
    echo 'module load intel' >> $bash_out
    echo 'module load curl' >> $bash_out
    echo 'module load java' >> $bash_out
    echo 'module load R/3.3.1' >> $bash_out
    #echo 'module load breseq/0.30' >> $bash_out
    echo 'module load breseq/0.27' >> $bash_out
    echo 'module load samtools/1.3.1' >> $bash_out
    echo '' >> $bash_out
    echo "breseq -j 8 -p -o ${OUT} -r ${REF} ${reads}" >> $bash_out

    OUT_annotated="${OUT}/output/evidence/annotated.gd"
    essentials_annotated="${OUT_essentials}/annotated.gd"
    echo "cp ${OUT_annotated} ${essentials_annotated}" >> $bash_out

    OUT_evidence="${OUT}/output/evidence/evidence.gd"
    essentials_evidence="${OUT_essentials}/evidence.gd"
    echo "cp ${OUT_evidence} ${essentials_evidence}" >> $bash_out

    OUT_output="${OUT}/output/output.gd"
    essentials_output="${OUT_essentials}/output.gd"
    echo "cp ${OUT_output} ${essentials_output}" >> $bash_out

    bam="${OUT}/data/reference.bam"
    coverage="${OUT}/data/coverage.txt"
    echo "samtools mpileup ${bam} > ${coverage}" >> $bash_out

    coverage_clean="${OUT_essentials}/coverage_clean.txt"
    mkdir -p "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/${DAY}/${line}"
    #echo "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/${DAY}/${line}"
    OUT_coverage="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_output_gbk_coverage/${DAY}/${line}/coverage_clean.txt"

    echo "python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/poly/getCoverage.py -c -i ${coverage} -o ${OUT_coverage}" >> $bash_out
    qsub $bash_out
  fi
done
