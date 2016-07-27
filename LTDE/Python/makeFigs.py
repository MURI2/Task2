from __future__ import division
import pandas as pd
import os, math, re
import numpy as np
import csv, collections
from itertools import chain
from scipy import stats
import  matplotlib.pyplot as plt
from operator import itemgetter
import scipy.stats as stats

mydir = os.path.expanduser("~/github/Task2/LTDE")


genus_color = {'Arthrobacter':'cyan', 'Janthinobacterium':'lightblue', 'Yersinia':'lightgreen', \
    'Rhodococcus':'tan', 'Bradyrhizobium':'pink', 'Pseudomonas':'darkgreen', 'Bacillus':'indigo'}

def plot_iRep(path = mydir):
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t' )
    strains = list(IN.Strain.unique())
    values = []
    genera = list(IN.Genus.unique())
    colors = []
    for genus in genera:
        values.append( IN[IN.Genus == genus]['iRep'].tolist())
        colors.append(genus_color[genus])
    fig = plt.figure(1, figsize=(9, 6))
    # Create an axes instance
    ax = fig.add_subplot(111)
    # Create the boxplot
    bp = ax.boxplot(values, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.3)
    ax.set_xticklabels(genera, fontsize = 8)
    anova = stats.f_oneway(*values)
    print "1-way ANOVA"
    print "F-value = " + str(anova[0])
    print "p-value = " + str(anova[1])
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_ylim(1, 2.7)
    ax.set_title('Genera')
    ax.set_ylabel('Average genome copy number')
    plt.plot([1, 1, 7, 7] , [2.4, 2.6, 2.6, 2.4], lw=1.5, c='k')
    plt.text(4, 2.6, "ns", ha='center', va='bottom', color='k', fontsize = 17)
    fig.savefig(mydir + '/figs/iRep.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def plotPopGenStats():
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    T_D = []
    params = ['Pi', 'W_theta', 'T_D']
    strains = list(IN.Strain.unique())
    genera = list(IN.Genus.unique())
    colors = []
    for genus in genera:
        colors.append(genus_color[genus])
    for param in params:
        param_values = []
        for strain in strains:
            param_values.append( IN[IN.Strain == strain][param].values)
        fig = plt.figure(1, figsize=(9, 6))
        # Create an axes instance
        ax = fig.add_subplot(111)
        ax.axhline(linewidth=2, color='darkgrey',ls='--')
        # Create the boxplot
        bp = ax.boxplot(param_values, patch_artist=True)

        #'maroon', 'blueviolet'
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.3)
        ax.set_xticklabels(strains)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        plt.legend()
        if param == 'T_D':
            ax.set_ylabel('Tajimas D')
            ax.set_ylim(-2, 4)
            plt.axhline(linewidth=2, color='darkgrey',ls='--')
        elif param == 'W_theta':
            ax.set_ylabel('Wattersons theta')
            ax.set_ylim(0, 2.5)
        elif param == 'Pi':
            ax.set_ylabel('pi')
            ax.set_ylim(0, 5)
        # Save the figure
        ax.set_title('Genera')
        ax.set_xticklabels(genera, fontsize = 8)
        fig.savefig(mydir + '/figs/' + param + '.png', bbox_inches='tight',  dpi = 600)
        plt.close()

def TDvsSlopePlot(evol = True):
    '''
    'evol' is the argument for if you only want the lines that fit the nonlinear model
    False means that all data is used.
    '''
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    x = IN.decay.values
    y = IN.T_D.values
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.decay, group.T_D, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10},  loc='upper right', frameon=False)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ylim([min(y)-0.1,max(y)+0.1])
    plt.xlim([min(x)-0.001,max(x)+0.001])
    plt.ylabel('Tajimas D', fontsize=20)
    plt.xlabel('Slope', fontsize=20)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/TD_vs_Slope' + evolTxt  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def SlopevsEvolPlot(evol = True):
    '''
    'evol' is the argument for if you only want the lines that fit the nonlinear model
    False means that all data is used.
    '''
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    x = IN.decay.values
    y = IN.evolvability.values
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.decay, group.evolvability, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10},  loc='upper right')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ylim([min(y)-0.000001,max(y)+0.000001])
    plt.xlim([min(x)-0.001,max(x)+ 0.001])
    plt.ylabel('Evolvability', fontsize=20)
    plt.xlabel('Slope', fontsize=20)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/Slope_vs_Evol' + evolTxt  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def EvolvsTDPlot(evol = True):
    '''
    'evol' is the argument for if you only want the lines that fit the nonlinear model
    False means that all data is used.
    '''
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    x = IN.evolvability.values
    y = IN.T_D.values
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.evolvability, group.T_D, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10},  loc='upper left')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.xlim([min(x)-0.000001,max(x)+0.000001])
    plt.ylim([min(y)-0.1,max(y)+ 0.1])
    plt.ylabel('Tajimas D', fontsize=20)
    plt.xlabel('Evolvability', fontsize=20)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/Slope_vs_Evol' + evolTxt  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def iRepvsSlopePlot(evol = True):
    '''
    'evol' is the argument for if you only want the lines that fit the nonlinear model
    False means that all data is used.
    '''
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    y = IN.decay.values
    x = IN.iRep.values
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.iRep, group.decay, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10})
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ylim([min(y)-0.001,max(y)+0.001])
    plt.xlim([min(x)-0.1,max(x)+ 0.1])
    plt.xlabel('Average genome copy number', fontsize=20)
    plt.ylabel('Slope', fontsize=20)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/iRep_vs_Slope' + evolTxt  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def iRepvsEvolPlot(evol = True):
    '''
    'evol' is the argument for if you only want the lines that fit the nonlinear model
    False means that all data is used.
    '''
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    x = IN.iRep.values
    y = IN.evolvability.values
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.iRep, group.evolvability, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10})
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ylim([min(y)-0.000001,max(y)+0.000001])
    plt.xlim([min(x)-0.1,max(x)+ 0.1])
    plt.ylabel('Evolvability', fontsize=20)
    plt.xlabel('Average genome copy number', fontsize=20)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    plt.savefig(mydir + '/figs/iRep_vs_Evol' + evolTxt  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


#plotPopGenStats()
#plot_iRep()
#TDvsSlopePlot(evol = True)
#TDvsSlopePlot(evol = False)
#SlopevsEvolPlot(evol = False)
#EvolvsTDPlot(evol = False)
#iRepvsSlopePlot(evol=True)
#iRepvsEvolPlot()
TDvsSlopePlot(evol = True)
