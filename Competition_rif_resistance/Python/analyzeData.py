from __future__ import division
import pandas as pd
import os, decimal
import numpy as np
import matplotlib.pyplot as plt

mydir = os.path.expanduser("~/GitHub/Task2/Competition/")


def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

def clean_data(strain = 'J'):
    IN_name = mydir + 'data/Competition_Ancestors_' + strain +'.csv'
    IN = open(IN_name, 'r')
    OUT_name = mydir + 'data/data_clean/Competition_Ancestors_' + strain + '.txt'
    OUT = open(OUT_name, 'w')
    print>> OUT, 'Exp_rep', 'Tech_rep', 'Anc0', 'Rif0', 'Anc1', 'Rif1'
    cell1 = 1
    cell2 = 1
    for line in IN.readlines():
        line = line.split(',')
        # get order of magnitudes for each of the two samples
        if '10E' in line[1]:
            cell1_str = find_between(line[1], '(*', ')')
            d1 = decimal.Decimal(cell1_str)
            cell1 = float(d1.to_eng_string())
            if strain == 'J':
                cell2_slice = 5
            else:
                cell2_slice = 4
            cell2_str = find_between(line[cell2_slice], '( ', ')')
            d2 = decimal.Decimal(cell2_str)
            cell2 =float(d2.to_eng_string())
        if '-' in line[0]:
            exp_tech = line[0].split('-')
            if line[3] == '':
                Anc0 = 'nan'
            else:
                Anc0 = int(line[3]) * cell1
            if line[2] == '':
                Rif0 = 'nan'
            else:
                Rif0 = int(line[2]) * cell1

            if strain == 'J':
                Anc1_slice = 7
                Rif1_slice = 6
            else:
                Anc1_slice = 6
                Rif1_slice = 5
            if line[Anc1_slice] == '':
                Anc1 = 'nan'
            else:
                #print int(line[Anc1_slice])
                Anc1 = int(line[Anc1_slice]) * cell2
            if line[Rif1_slice] == '':
                Rif1 = 'nan'
            else:
                Rif1 = int(line[Rif1_slice]) * cell2
            #Anc1 = int(line[6]) * cell2
            #Rif1 = int(line[7]) * cell2
            # we only want the small CPU and CPU on rif
            print>> OUT, exp_tech[0], exp_tech[1], Anc0, Rif0, Anc1, Rif1
    OUT.close()


def r(row):
    #print np.log(row['Rif1'] / row['Rif0']), np.log(row['Anc1'] / row['Anc0'])
    return np.log(row['Rif1'] / row['Rif0']) - np.log(row['Anc1'] / row['Anc0'])



def plotr():
    strains = ['B', 'D', 'J', 'P']
    strains_names = ['Bacillus', 'Deinococcus', 'Janthino', 'Pseudomonas']
    strain_lists = []
    for strain in strains:
        IN_file = mydir + 'data/data_clean/Competition_Ancestors_' + strain + '.txt'
        IN = pd.read_csv(IN_file, sep = ' ')
        IN.dropna()
        IN = IN.dropna()
        IN['r'] = IN.apply(r, axis=1)
        #IN_r = IN[['Exp_rep','r']]
        #print IN_r
        IN_group = IN.groupby(['Exp_rep']).mean()
        print IN['r']
        strain_lists.append(IN_group['r'].values)
    # Create a figure instance
    fig = plt.figure(1, figsize=(9, 6))

    # Create an axes instance
    ax = fig.add_subplot(111)
    # Create the boxplot
    bp = ax.boxplot(strain_lists)
    ax.axhline(y=0, linewidth=2, color='grey', linestyle = '--')
    ax.set_ylabel('Selection rate', fontsize = 18)
    ax.set_title('Selection rate between rifampicin resistant strain and ancestor \n' + r'$r_{Rif} - r_{Ancestor}$')
    ## Custom x-axis labels
    ax.set_xticklabels(strains_names)
    fig.savefig(mydir + 'figs/fig1.png', bbox_inches='tight', dpi=600)
    plt.close()


strains = ['B', 'D', 'J', 'P']
for strain in strains:
    clean_data(strain = strain)
plotr()
