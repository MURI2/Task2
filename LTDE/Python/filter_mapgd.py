from __future__ import division
import pandas as pd
import os
import numpy as np
import scipy.stats as st

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'

IN = pd.read_csv(mydir + 'mapgd/annotate/KBS0727_merged_annotate.pol', delimiter = '\t', \
    header = None)

def pandasNames(df):
    names = ['Scaffold', 'Pos', 'Major', 'Minor', 'Coverage', 'Error', \
        'Gene', 'AA-change']
    samples = len(df.columns) - len(names)
    for x in range(samples):
        name = 'Sample_' + str(x+1)
        names.append(name)
    df.columns = names
    print df.iloc[:,8:]
    #df.iloc[:,8:] = pd.to_numeric(df.iloc[:,8:], errors='coerce')
    df.iloc[:,8:] = df.iloc[:,8:].convert_objects(convert_numeric=True)
    #df.iloc[:,8:] = df.iloc[:,8:].astype(float)
    return (df, samples)

names_IN = pandasNames(IN)
IN = names_IN[0]
IN_samples = names_IN[1]


ones_and_zeroes = IN[IN.apply(lambda x: ((x.iloc[8:]==float(1)).sum() !=0 \
        or (x.iloc[8:]==float(0)).sum() !=0), axis=1) ]

#ones_and_zeroes = IN[IN.apply(lambda x: ((x.iloc[8:]==float(1)).sum() !=0 \
#        or (x.iloc[8:]==float(0)).sum() !=0) and  x.iloc[6] != 'NC', axis=1) ]

print ones_and_zeroes
#print ones_and_zeroes.iloc[:,8:]
