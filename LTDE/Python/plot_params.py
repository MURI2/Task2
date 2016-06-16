from __future__ import division
import pandas as pd
import os, math
import numpy as np
import csv, collections
from itertools import chain
from scipy import stats
import  matplotlib.pyplot as plt
from operator import itemgetter

mydir = os.path.expanduser("~/github/Task2/LTDE")

def plotPopGenStats():
    IN = pd.read_csv(mydir + '/data/mapgd/final/PopGenStats.txt', sep = ' ', header = None)
    IN.columns = ['Strain', 'Rep', 'Pi', 'Pi_Var', 'W', 'W_Var', 'T_D']
    #strains = []
    T_D = []
    params = ['Pi', 'W', 'T_D']
    for param in params:
        strains = IN.Strain.unique()
        param_values = []
        for strain in strains:
            #IN[['Pi']]surveys_df[surveys_df.year == 2002]
            param_values.append( IN[IN.Strain == strain][param].values)

        fig = plt.figure(1, figsize=(9, 6))
        # Create an axes instance
        ax = fig.add_subplot(111)
        ax.axhline(linewidth=2, color='darkgrey',ls='--')
        # Create the boxplot
        bp = ax.boxplot(param_values, patch_artist=True)
        colors = ['cyan', 'lightblue', 'lightgreen', 'tan', 'pink', 'darkgreen', \
            'maroon', 'blueviolet', 'indigo']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.3)

        ax.set_xticklabels(strains)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        if param == 'T_D':
            ax.set_ylabel('Tajimas D')
            ax.set_ylim(0, 4)
        elif param == 'W':
            ax.set_ylabel('Wattersons theta')
            ax.set_ylim(0, 4)
        elif param == 'pi':
            ax.set_ylabel('pi')
            #plt.ylabel('ppiiii')
            ax.set_ylim(0, 4)
        # Save the figure
        ax.set_title('Strains')
        fig.savefig(mydir + '/figs/' + param + '.png', bbox_inches='tight',  dpi = 600)
        plt.close()


def TDvsEvolData():
    IN_TD = (mydir + '/data/mapgd/final/PopGenStats.txt')
    IN_Evol = (mydir + '/data/perRepDeathCurveTraits.txt')
    strains = []
    data = []
    test_dict = {}
    with open(IN_TD) as f:
        my_lines = f.readlines()
        for x in my_lines:
            x = x.strip().split(' ')
            if x[0] == 'KBS0721' or x[0] == 'KBS0710':
                continue
            strains.append(x[0])
            #data_x = [ float(y) for y in x[6] ]
            data_x = float(x[6])
            if x[0] in test_dict:
                pass
            else:
                test_dict[x[0]] = {}

            test_dict[x[0]][int(x[1])] = [data_x]
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
                    if x[1] == '10':
                        test_dict[x[0]][4].extend([x[2], x[3]])
                    elif x[1] == '11':
                        test_dict[x[0]][5].extend([x[2], x[3]])
                    #elif x[1] == '1':
                    #    test_dict[x[0]][].extend([x[2], x[3]])
                    elif x[1] == '3':
                        test_dict[x[0]][2].extend([x[2], x[3]])
                    elif x[1] == '4':
                        test_dict[x[0]][3].extend([x[2], x[3]])
                    else:
                        continue
                else:
                    test_dict['KBS0711'][int(x[1]) + 5].extend([x[2], x[3]])
    return test_dict



def get_list(d):
    return_list  =[]
    for key, value in d.iteritems():
        for nested_key, nested_value in value.iteritems():
            nested_tuple = (key, nested_value)
            return_list.append( nested_tuple)
            #return_list.append( nested_value)
    return return_list

def TDvsEvolPlot():
    data = TDvsEvolData()
    nested_list =get_list(data)
    x = np.asarray([z[1][0] for z in nested_list ])
    y = np.asarray([z[1][1] for z in nested_list ])
    evol = [abs(z[1][2]) for z in nested_list ]
    names = np.asarray([z[0] for z in nested_list ])
    names_count = collections.Counter(names)
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    colors = ['cyan', 'lightblue', 'lightgreen', 'tan', 'pink', 'darkgreen', \
        'maroon', 'blueviolet', 'indigo']
    colors_count = []
    for enum, value in enumerate(names_count_sorted):
        color = colors[enum]
        colors_count.extend([color for i in range(value[1])])
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y, c= colors_count, alpha = 0.8, s = 48,edgecolor='black')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.axhline(linewidth=2, color='darkgrey',ls='--')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.xlim([0.5,3.5])
    plt.ylim([-0.005,0.020])
    plt.xlabel('Tajimas D', fontsize=20)
    plt.ylabel('Slope', fontsize=20)
    fig.savefig(mydir + '/figs/test.png',  bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

plotPopGenStats()

#TDvsEvolData()
#TDvsEvolPlot()
