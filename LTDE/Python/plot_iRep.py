from __future__ import division
import pandas as pd
import os, math, re
import numpy as np
import csv, collections
from itertools import chain
from scipy import stats
import  matplotlib.pyplot as plt

mydir = os.path.expanduser("~/github/Task2/LTDE")

def clean_iRep(path = mydir):
    strains = ['KBS0703', 'KBS0710', 'KBS0711', 'KBS0713', \
        'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802', 'KBS0812']
    lines_to_keep = [6, 10, 14, 18]
    to_pandas = []
    for strain in strains:
        strainPath = path + '/data/iRep/' + strain + '/'
        for content in os.listdir(strainPath):
            if content.endswith('.tsv'):
                # iRep, r2, coverage, % windows passing filer.
                name = re.split(r'[._]+', content)
                content_list = []
                content_list.append(name[0])
                content_list.append(name[1])
                for count, line in enumerate(open(strainPath + content)):
                    split_line = re.split(r'[#\t\n]+', line)
                    # remove empty items
                    clean_line = filter(None, split_line)
                    if len(clean_line) == 0:
                        continue
                    if count in lines_to_keep:
                        content_list.append(float(clean_line[-1]))
                to_pandas.append( content_list)
    df = pd.DataFrame(to_pandas)
    df.columns  =['Strain', 'Replicate', 'iRep', 'R2', 'Coverage', 'PercentPass']
    df_path = path + '/data/iRep/iRepFinal.txt'
    df.to_csv(df_path, sep='\t', index = False)

def plot_iRep(path = mydir):
    IN = pd.read_csv(mydir + '/data/iRep/iRepFinal.txt', sep ='\t' )
    strains = list(IN.Strain.unique())
    values = []
    for strain in strains:
        values.append( IN[IN.Strain == strain]['iRep'].values)
    fig = plt.figure(1, figsize=(9, 6))
    # Create an axes instance
    ax = fig.add_subplot(111)
    bp = ax.boxplot(values)
    ax.set_xticklabels(strains, fontsize = 8)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_ylim(1, 2.2)
    ax.set_title('Strains')
    ax.set_ylabel('Index of replication')
    fig.savefig(mydir + '/figs/iRep.png', bbox_inches='tight',  dpi = 600)
    plt.close()

plot_iRep()
