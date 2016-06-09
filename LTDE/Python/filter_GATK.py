from __future__ import division
import pandas as pd
import os, math
import numpy as np
import scipy.stats as st
import csv, re
mydir = os.path.expanduser("~/github/Task2/LTDE")

taxa = ['KBS0703', 'KBS0705', 'KBS0706', 'KBS0710', 'KBS0711', 'KBS0713', \
    'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802']


def getSNPTable(strains):
    for strain in strains:
        dfs = []
        samples = []
        path = mydir + '/data/GATK/annotate/' + strain
        for i in os.listdir(path):
            if i.endswith(".txt"):
                sample_name = i.split('.')[0]
                sample = re.split(r'[-_]+', sample_name)[2]
                samples.append(sample)
                names = ['Scaffold', 'Pos', 'Ref', 'Alt', 'Gene', \
                    'Coding', 'Qual', 'AC']
                IN = pd.read_csv(path + '/' + i, delimiter = ' ', header = None)
                IN.columns = names
                # remove original header and AC column
                IN = IN.iloc[1:,:-1]
                dfs.append(IN)
            else:
                continue
        merged_df = reduce(lambda left,right: pd.merge(left,right,on=['Scaffold', 'Pos','Ref','Gene'], how='outer'), dfs)
        toRename = ['Alt', 'Coding', 'Qual']
        new_columns = ['Scaffold', 'Pos', 'Ref', 'Gene']
        to_move = ''
        for x, y in enumerate(samples):
            for z in toRename:
                out = z + '_' + y
                if x ==0 and z == 'Alt':
                    new_columns.insert(3,out)
                    to_move = out
                else:
                    new_columns.append(out)
        merged_df.columns = new_columns
        to_move_values = merged_df.loc[:,to_move]
        merged_df.drop(to_move, axis=1, inplace=True)
        merged_df.insert(4, to_move, to_move_values)
        merged_df.to_csv(mydir + '/data/GATK/merged/all/' + strain +'_GATK.txt' ,sep='\t', \
            index = False)

def countNANs(row):
    row_list =  np.asarray(list(row))
    print row.value_counts()
    print row

def filterSNPTable(strains, coding_to_csv=False):
    for strain in strains:
        path = mydir + '/data/GATK/merged/all/' + strain + '_GATK.txt'
        IN = pd.read_csv(path, delimiter = '\t', header = 'infer')
        IN_coding = IN[IN['Gene'] != 'NC']
        if coding_to_csv == True:
            IN_coding.to_csv(mydir + '/data/GATK/merged/coding/' + strain +'_GATK_C.txt' ,sep='\t', \
                index = False)
        NAN_counts = IN_coding.isnull().sum(axis=1)
        toKeep = NAN_counts[NAN_counts != 0].index
        toKeep_df = IN_coding.ix[toKeep]
        toKeep_df.to_csv(mydir + '/data/GATK/merged/potential_coding_mutations/' + strain +'_GATK_C_PCM.txt' ,sep='\t', \
            index = False)


#getSNPTable(taxa)
#filterSNPTable(taxa)
