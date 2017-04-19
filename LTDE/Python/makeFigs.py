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


genus_color = {'Arthrobacter':'cyan', 'Janthinobacterium':'lightblue', \
        'Yersinia':'lightgreen', 'Rhodococcus':'tan', 'Bradyrhizobium':'pink', \
        'Pseudomonas':'darkgreen', 'Bacillus':'indigo', 'Variovorax': 'darkred',\
        'Burkholderia': 'orange', 'Curtobacterium': 'red'}

def plot_iRep(evol = True, path = mydir):
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t' )
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    else:
        IN = IN.loc[IN['evolvability'] >= float(0)]
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
    bp = ax.boxplot(values, patch_artist=True, bootstrap = 1000)
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
    ax.set_ylim(1.15, 1.7)
    #ax.set_title('Genera')
    ax.set_ylabel('Average genome copy number')
    #plt.plot([1, 1, 7, 7] , [1.6, 1.65, 1.65, 1.6], lw=1.5, c='k')
    #plt.text(4, 1.66, "****", ha='center', va='bottom', color='k', fontsize = 17)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/iRep/iRep' + evolTxt +'.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def plotPopGenStats(evol = True):
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    else:
        IN = IN.loc[IN['evolvability'] >= float(0)]
    T_D = []
    params = ['Pi', 'W_theta', 'T_D']
    #params = ['T_D']
    strains = list(IN.Strain.unique())
    genera = list(IN.Genus.unique())
    colors = []
    values = []
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
        bp = ax.boxplot(param_values, patch_artist=True, bootstrap = 1000)
        #'maroon', 'blueviolet'
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.3)
        ax.set_xticklabels(strains)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        plt.legend()
        if param == 'T_D':
            ax.set_ylabel(r'$D_{T}$', fontsize = 26)
            ax.set_ylim(0, 4)
            plt.axhline(linewidth=2, color='darkgrey',ls='--')
            #plt.plot([1, 1, 7, 7] , [3.5, 3.7, 3.7, 3.5], lw=1.5, c='k')
            #plt.text(4, 3.73, "****", ha='center', va='bottom', color='k', fontsize = 17)
            anova = stats.f_oneway(*param_values)
            print "1-way ANOVA"
            print "F-value = " + str(anova[0])
            print "p-value = " + str(anova[1])
        elif param == 'W_theta':
            #ax.set_ylabel('Wattersons theta')
            ax.set_ylabel(r'$\hat{\theta}_{W}$')
            ax.set_ylim(0, 2.5)
        elif param == 'Pi':
            #ax.set_ylabel('pi')
            ax.set_ylabel(r'$\pi$')
            ax.set_ylim(0, 5)
        # Save the figure
        ax.set_title('Genera')
        ax.set_xticklabels(genera, fontsize = 8)
        if evol == True:
            evolTxt = '_evol'
        else:
            evolTxt = ''
        fig.savefig(mydir + '/figs/PopGen/' + param + evolTxt+ '.png', bbox_inches='tight',  dpi = 600)
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
    else:
        IN = IN.loc[IN['evolvability'] >= float(0)]
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
    ax.legend(numpoints=1, prop={'size':10},  loc='upper left', frameon=False)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    #plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ylim([min(y)-0.1,max(y)+0.1])
    plt.xlim([min(x)-0.001,max(x)+0.001])
    plt.ylabel(r'$D_{T}$' + ' (lack of rare alleles)', fontsize=20)
    plt.xlabel('Growth rate', fontsize=20)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/TD_vs_Slope' + evolTxt  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def iRepvsTDPlot(evol = True):
    '''
    'evol' is the argument for if you only want the lines that fit the nonlinear model
    False means that all data is used.
    '''
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    else:
        IN = IN.loc[IN['evolvability'] >= float(0)]
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    x = IN.iRep.values
    y = IN.T_D.values
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.iRep, group.T_D, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10},  loc='upper left', frameon=False)
    ax.text(1.2, 2.68, r'$r^{2}=0.751$', fontsize=14)
    ax.text(1.2, 2.5 , r'$p < 0.05$', fontsize=14)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "intecept = "  + str(intercept)
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
    plt.xlim([min(x)-0.05,max(x)+ 0.05])

    plt.ylabel(r'$D_{T}$' + ' (lack of rare alleles)', fontsize=20)
    plt.xlabel('Average genome copy number', fontsize=20)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/iRep_vs_TD' + evolTxt  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def iRepvsPi(evol = True):
    '''
    'evol' is the argument for if you only want the lines that fit the nonlinear model
    False means that all data is used.
    '''
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
        IN = IN.loc[IN['iRep'] < float(7)]
    else:
        IN = IN.loc[IN['evolvability']  >= float(0)]
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    x = IN.iRep.values
    y = IN.Pi.values
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.iRep, group.Pi, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10},  loc='upper left', frameon=False)
    ax.text(1.2, 2.68, r'$r^{2}=0.751$', fontsize=14)
    ax.text(1.2, 2.5 , r'$p < 0.05$', fontsize=14)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "intecept = "  + str(intercept)
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
    plt.xlim([min(x)-0.3,max(x)+ 0.05])

    plt.ylabel(r'$\pi$', fontsize=20)
    plt.xlabel('Average genome copy number', fontsize=20)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/iRep_vs_Pi' + evolTxt  + '.png', \
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
    ax.legend(numpoints=1, prop={'size':10},  loc='upper right', frameon=False)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    #plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ylim([min(y)-0.000001,max(y)+0.000001])
    plt.xlim([min(x)-0.001,max(x)+ 0.001])
    plt.ylabel('Change in growth rate', fontsize=20)
    plt.xlabel('Decay', fontsize=20)
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
    pd.set_option('display.precision',20)
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    x = IN.evolvability.values * 2
    y = IN.T_D.values
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.evolvability * 2, group.T_D, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10},  loc='upper left', frameon=False)
    ax.text(0.000006, 2.58, r'$r^{2}=0.635$', fontsize=14)
    ax.text(0.000006, 2.4 , r'$p < 0.001$', fontsize=14)
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
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.xlim([min(x)-0.000001,max(x)+0.000001])
    plt.ylim([min(y)-0.1,max(y)+ 0.1])
    plt.ylabel('Population structure ' + r'$(D_{T})$' , fontsize=24)
    plt.xlabel('Change in growth rate', fontsize=24)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.savefig(mydir + '/figs/Evol_vs_TD' + evolTxt  + '.png', \
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
        IN = IN.loc[IN['iRep'] < float(7)]
    else:
        IN = IN.loc[IN['evolvability'] >= float(0)]
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
    ax.legend(numpoints=1, prop={'size':10}, frameon=False, loc='lower left')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    #plt.plot(x, predict_y, 'k-')
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
        IN = IN.loc[IN['iRep'] < float(7)]
    x = IN.iRep.values
    y = IN.evolvability.values *2
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.iRep, group.evolvability *2, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10}, frameon=False, loc='upper left')
    ax.text(1.15, 0.0000165, r'$r^{2}=0.692$', fontsize=14)
    ax.text(1.15, 0.0000155, r'$p \:\ll  0.0001$', fontsize=14)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "intecept = "  + str(intercept)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylim([min(y)-0.000001,max(y)+0.000001])
    plt.xlim([min(x)-0.1,max(x)+ 0.1])
    plt.ylabel('Change in growth rate', fontsize=24)
    plt.xlabel('Average genome copy number \n (proxy for birth rate)', fontsize=24)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    plt.savefig(mydir + '/figs/iRep_vs_Evol' + evolTxt  + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def PivsEvolPlot(evol = True):
    '''
    'evol' is the argument for if you only want the lines that fit the nonlinear model
    False means that all data is used.
    '''
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
        IN = IN.loc[IN['iRep'] < float(7)]
    x = IN.evolvability.values *2
    y = IN.Pi_prob.values
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    groups = IN.groupby('Genus')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(group.evolvability *2, group.Pi_prob, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=name, c = genus_color[name])
    ax.legend(numpoints=1, prop={'size':10}, frameon=False, loc='upper left')
    ax.text(min(x), min(y) + 0.0000016, r'$r^{2}=0.527$', fontsize=14)
    ax.text(min(x), min(y) + 0.00000145, r'$p \ll \, 0.05$', fontsize=14)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "intecept = "  + str(intercept)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    plt.ylabel('Nucleotide diversity, ' + r'$(\pi)$', fontsize=22, labelpad=20)
    plt.xlabel('Change in growth rate', fontsize=22)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    plt.savefig(mydir + '/figs/PivsEvolPlot' + evolTxt  + '.png', \
          pad_inches = 0.4, dpi = 600)
    plt.close()


def SlopevsEvolvsPiPlot(evol = True):
    IN = pd.read_csv(mydir + '/data/Final/MAPGD_Evol_iRep.txt', sep ='\t')
    pd.set_option('display.precision', 20)
    if evol == True:
        IN = IN.loc[IN['evolvability'] > float(0)]
        IN = IN.loc[IN['iRep'] < float(7)]

    fig = plt.figure()
    count = 0
    for i in range(2):
        ax = fig.add_subplot(2, 1, i+1)
        y = IN.T_D.values
        if i == 0:
            x = IN.decay.values
        elif i == 1:
            x = IN.evolvability.values *2
        names_count = collections.Counter(IN.Strain.tolist())
        names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
        groups = IN.groupby('Genus')

        ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
        for name, group in groups:
            if i == 0:
                ax.plot(group.decay, group.T_D, marker='o', alpha = 0.8, \
                    linestyle='', ms=12, label=name, c = genus_color[name])
            elif i == 1:
                ax.plot(group.evolvability *2, group.T_D, marker='o', alpha = 0.8, \
                    linestyle='', ms=12, label=name, c = genus_color[name])
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        print "intecept = "  + str(intercept)
        print "slope = " + str(slope)
        print "r2 = " + str(r_value**2)
        print "p = " + str(p_value)
        predict_y = intercept + slope * x
        pred_error = y - predict_y
        degrees_of_freedom = len(x) - 2
        residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

        if i == 1:
            plt.plot(x , predict_y, 'k-')
            ax.text(min(x), max(y) - 0.3, r'$r^{2}=0.641$', fontsize=14)
            ax.text(min(x), max(y) - 0.65 , r'$p \:\ll   0.001$', fontsize=14)
            ax.set_xlabel('Change in growth rate', fontsize = 16)
        elif i == 0:
            ax.legend(numpoints=1, prop={'size':10}, frameon=False, loc='upper left')
            ax.text(max(x) - 0.0004, max(y) - 0.3, r'$r^{2}=0.049$', fontsize=14)
            ax.text(max(x) - 0.0004, max(y) - 0.65, r'$p \:\nless  0.05$', fontsize=14)
            ax.set_xlabel('Growth rate', fontsize = 16)
    if evol == True:
        evolTxt = '_evol'
    else:
        evolTxt = ''
    fig.subplots_adjust(hspace=0.24)
    fig.text(0.05, 0.5, r'$D_{T}$', ha='center', va='center', rotation='horizontal', fontsize=22)
    plt.savefig(mydir + '/figs/SlopevsEvolvsPiPlot' + evolTxt  + '.png', \
          pad_inches = 0.4, dpi = 600)
    plt.close()



#iRepvsEvolPlot()
#EvolvsTDPlot()
#iRepvsPi()
#iRepvsSlopePlot()
#TDvsSlopePlot()
#SlopevsEvolvsPiPlot()
#PivsEvolPlot()

Fig3()
