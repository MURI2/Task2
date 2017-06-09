from __future__ import division
import cleanData as cd
import pandas as pd
import os
import scipy.stats as stats
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.cov_struct import (Exchangeable,
    Independence,Autoregressive)
from statsmodels.genmod.families import Poisson

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/data/")

#sample_by_gene_matrix('D100', 'P')
#run_everything('D100', split = True, get_variants = True, merge_variants = True,\
#    unique_mutations = True, multiple_gene_hits = True, sample_by_gene_matrix = True, \
#    variant_type = 'DEL')
def cleanGBK(strans):
    for strain in strains:
        print strain
        cd.cleanGBK(strain)

strains = ['B', 'C', 'D', 'F', 'J', 'P']
#cleanGBK(strains)

def poisRegress():
    IN_file_mut = mydir + 'gene_by_sample/D/D100/sample_by_gene.txt'
    IN_file_genes = mydir + 'reference_assemblies_task2_table/D.txt'
    IN_mut = pd.read_csv(IN_file_mut, sep = '\t')
    IN_genes = pd.read_csv(IN_file_genes, sep = ' ')
    print IN_genes

poisRegress()
