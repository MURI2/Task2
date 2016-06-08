from __future__ import division
import pandas as pd
import os, math
import numpy as np
import csv, collections
from itertools import chain
from scipy import stats
import  matplotlib.pyplot as plt

mydir = os.path.expanduser("~/github/Task2/LTDE")

def plotPopGenStats():
    IN = pd.read_csv(mydir + '/data/mapgd/final/PopGenStats.txt', sep = ' ', header = None)
    IN.columns = ['Strain', 'Rep', 'Pi', 'Pi_Var', 'W', 'W_Var', 'T_D']
    #strains = []
    T_D = []
    params = ['Pi', 'W', 'T_D']
    strains = IN.Strain.unique()
    for param in params:
        param_values = []
        for strain in strains:
            #IN[['Pi']]surveys_df[surveys_df.year == 2002]
            param_values.append( IN[IN.Strain == strain][param].values)
        fig = plt.figure(1, figsize=(9, 6))
        # Create an axes instance
        ax = fig.add_subplot(111)
        ax.axhline(linewidth=2, color='darkgrey',ls='--')
        # Create the boxplot
        bp = ax.boxplot(param_values)
        ax.set_xticklabels(strains)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        if param == 'T_D':
            ax.set_ylabel('Tajimas D')
            ax.set_ylim(-3.5, 3.5)
        elif param == 'W':
            ax.set_ylabel('Wattersons theta')
            print param_values
            ax.set_ylim(0, 3.5)
        elif param == 'pi':
            ax.set_ylabel('nucleotide diversity (pi)')
            ax.set_ylim(0, 3.5)
        # Save the figure
        fig.savefig(mydir + '/figs/' + param + '.png', bbox_inches='tight')
        plt.close()


#plotTD()

strain_list = ['KBS0722', 'KBS0724', 'KBS0727', 'KBS0715', 'KBS0703', 'KBS0711', \
            'KBS0713', 'KBS0802']

def TDvsEvol():
    IN_TD = (mydir + '/data/mapgd/final/mapgd_TD.txt')
    IN_Evol = (mydir + '/data/perRepDeathCurveTraits.txt')
    strains = []
    data = []
    test_dict = {}
    with open(IN_TD) as f:
        my_lines = f.readlines()
        for x in my_lines:
            #None
            x = x.strip().split(',')
            if x[0] == 'KBS0721' or x[0] == 'KBS0710':
                continue
            strains.append(x[0])
            data_x = [ float(y) for y in x[1:] ]
            test_dict[x[0]] = {}
            for y, z in enumerate(data_x):
                y = y +1
                test_dict[x[0]][y] = [z]
            #print np.mean(data_x)
            data.append(data_x)
    with open(IN_Evol) as f:
        my_lines = f.readlines()
        for x in my_lines:
            x = x.strip().split(' ')
            if x[0] in strain_list or x[0] == 'KBS0711W':
                x[2] = float(x[2])
                x[3] = float(x[3])

                if x[0] != 'KBS0711' and x[0] != 'KBS0711W':
                    if int(x[1]) in test_dict[x[0]].keys():
                        #print test_dict[x[0]][int(x[1])]
                        test_dict[x[0]][int(x[1])].extend([x[2], x[3]])
                        #print test_dict[x[0]][int(x[1])]
                elif x[0] == 'KBS0711':
                    if x[1] == '11':
                        test_dict[x[0]][4].extend([x[2], x[3]])
                    elif x[1] == '1':
                        test_dict[x[0]][int(x[1])].extend([x[2], x[3]])
                    elif x[1] == '3':
                        test_dict[x[0]][int(2)].extend([x[2], x[3]])
                    elif x[1] == '4':
                        test_dict[x[0]][int(3)].extend([x[2], x[3]])
                else:
                    test_dict['KBS0711'][int(x[1]) + 4].extend([x[2], x[3]])

    return test_dict
    #print data



plotPopGenStats()

#test = TDvsEvol()
#print test.values()[1]

#def get_list(d):
#    return_list  =[]
#    for key, value in d.iteritems():
#        for nested_key, nested_value in value.iteritems():
#            return_list.append( nested_value)
#    return return_list

#nested_list =get_list(test)
#TD = [x[0] for x in nested_list]
#slope = [x[1] for x in nested_list]
#evol = [abs(x[2]) for x in nested_list]

#print nested_list
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(TD, slope, color='blue',s=5, edgecolor='none')
#slope, intercept, r_value, p_value, std_err = stats.linregress(TD,slope)
#print slope, r_value, p_value
#plt.xlim([0,100])
#plt.xlabel('Tajimas D', fontsize=20)
#plt.ylabel('Slope', fontsize=20)
#fig.savefig(mydir + '/figs/test.png',  bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
