from __future__ import division
import pandas as pd
import os
import numpy as np
import scipy.stats as st

mydir = os.path.expanduser("~/github/Task2/LTDE")

taxa = ['KBS0703', 'KBS0710', 'KBS0711', 'KBS0713', 'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802']

def pandasNames(df):
    names = ['Scaffold', 'Pos', 'Major', 'Minor', 'Coverage', 'Error', \
        'Gene', 'AA-change']
    samples = len(df.columns) - len(names)
    for x in range(samples):
        name = 'Sample_' + str(x+1)
        names.append(name)
    df.columns = names
    #df.iloc[:,8:] = pd.to_numeric(df.iloc[:,8:], errors='coerce')
    df.iloc[:,8:] = df.iloc[:,8:].convert_objects(convert_numeric=True)
    #df.iloc[:,8:] = df.iloc[:,8:].astype(float)
    return (df, samples)

def func(x):
    count = 0
    for y in x:
        if y != float(0) or y != float(1):
            count += 1

def mapgdPI(row):
    print row

def getPiVar(pi, n):
    '''
    returns the variance for pi
    got here: http://genapps.uchicago.edu/slider/index.html
    Also in Walsh & Lynch 201?
    '''
    if pi == 0:
        return 0
    else:
        x1 = ((n+1) / (3* (n-1))) * pi
        x2 = ((( (n**2) + n +3 ) * 2 ) / (9*n * (n-1)) ) * (pi**2)
        return x1 + x2



def getPolyTable(strains):
    polyDict = {}
    for taxon in taxa:
        IN = pd.read_csv(mydir + '/data/mapgd/annotate/' + taxon + '_merged_annotate.pol', delimiter = '\t', \
            header = None)
        names_IN = pandasNames(IN)
        IN = names_IN[0]
        IN_samples = names_IN[1]
        # count polymorphisms in coding vs noncoding

        # non coding and coding subsets
        NC_count = IN[IN.apply(lambda x: (x.iloc[6] == 'NC'), axis=1 )]
        C_count = IN[IN.apply(lambda x: (x.iloc[6] != 'NC'), axis=1 )]
        C_S_count = IN[IN.apply(lambda x:  (x.iloc[7] == 'S'), axis=1 )]
        samples_NC = []
        samples_C = []
        for x in range(1, IN_samples):
            sample_column = 'Sample_' + str(x)
            NC_x =  sum(sum([(NC_count[sample_column] > (float(0))) & (NC_count[sample_column]< (float(1)))]))
            #NC_x = (NC_count[sample_column]  > (float(0))) and  (NC_count[sample_column]  < (float(1)))
            C_x =  sum(sum([(C_count[sample_column] > (float(0))) & (C_count[sample_column]< (float(1)))]))
            #C_x = (C_count[sample_column]  != (float(0)  or float(1))).count()
            samples_NC.append(NC_x)
            samples_C.append(C_x)
        NC_mean = np.mean(samples_NC)
        NC_std = np.std(samples_NC)
        C_mean = np.mean(samples_C)
        C_std = np.std(samples_C)
        polyDict[taxon] = [NC_mean, NC_std, C_mean, C_std]
        #ones_and_zeroes = IN[IN.apply(lambda x: ((x.iloc[8:]==float(1)).sum() !=0 \
        #        or (x.iloc[8:]==float(0)).sum() !=0) and  x.iloc[6] != 'NC', axis=1) ]
        #ones_and_zeroes_path = mydir + '/data/mapgd/final/coding/' + taxon + '_coding.txt'
        #ones_and_zeroes.to_csv(ones_and_zeroes_path,sep='\t')
    df = pd.DataFrame(polyDict)
    df = df.transpose()
    df.columns = ['NC_mean', 'NC_std', 'C_mean', 'C_std']
    latex_path = mydir + '/tables/test.tex'
    df.to_latex(latex_path)

    #print df


def getPi(strains):
    piDict = {}
    for strain in strains:
        IN = pd.read_csv(mydir + '/data/mapgd/annotate/' + strain + '_merged_annotate.pol', delimiter = '\t', \
            header = None)
        names_IN = pandasNames(IN)
        IN = names_IN[0]
        IN_samples = names_IN[1]
        C_S_count = IN[IN.apply(lambda x:  (x.iloc[7] == 'S'), axis=1 )]
        pi = []
        piVar = []
        for x in range(1, IN_samples + 1):
            sample_column = 'Sample_' + str(x)
            pi_subset_names = ['Coverage', sample_column]
            pi_subset = C_S_count[pi_subset_names]
            pi_subset_tuples = [tuple(x) for x in pi_subset.values]
            pi_x = []
            N_x = []
            for value in pi_subset_tuples:
                if len(value) == 0 or value[1] == float(1) or value[1] == float(0):
                    continue
                p = value[1]
                q = 1- p
                N = value[0]
                pi_value = (N / (N-1) ) * p * q * 2
                pi_x.append(pi_value)
                N_x.append(N)
            # sum over pi
            pi_x_sum = sum(pi_x)
            # calculate mean N for the variance of pi
            N_x_mean = np.mean(N_x)
            pi_x_var = getPiVar(pi_x_sum, N_x_mean)
            pi.append(pi_x_sum)
            piVar.append(pi_x_var)
        pi_mean = np.mean(pi)
        pi_var = np.mean(piVar)
        piDict[strain] = [pi_mean, pi_var]

    df = pd.DataFrame(piDict, index=None)
    df = df.transpose()
    df.reset_index(level=0, inplace=True)
    df.columns = ['strain', 'pi_mean', 'pi_var']
    df_path = mydir + '/data/mapgd/final/mapgd_pi.txt'
    df.to_csv(df_path,sep='\t')





#getPolyTable(taxa)
getPi(taxa)
