from __future__ import division
import os, math, numbers, itertools
import pandas as pd
import numpy as np
from string import maketrans
from collections import Counter
import  matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/")

strain_colors = {'Janthinobacterium':'cyan', 'Caulobacter':'lightblue', \
         'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'indigo',
         'Pedobacter': 'darkred'}

species_dict = {'B': 'Bacillus', 'C':'Caulobacter', 'D':'Deinococcus', \
    'F': 'Pedobacter', 'J': 'Janthinobacterium', 'P':'Pseudomonas'}

treatment_dict = {'0': '1-day', '1': '10-day', '2': '100-day'}

def CV_KDE(oneD_array):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(kernel='exponential'),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple


def AFS():
    strains = ['B', 'C', 'D', 'F', 'J', 'P']
    #strains = ['B']
    colors = {'1':'cyan', '2':'lightblue', \
             '3': 'red', '4':'darkgreen', '5':'indigo'}
    treatments = ['0', '1', '2']
    for strain in strains:
        path = mydir + 'data/breseq_output_gbk_essentials_split_clean_merged_unique/D100/Strain_' + strain +'.txt'
        if os.path.exists(path) != True:
            continue
        IN = pd.read_csv(path, sep = '\t', header = 'infer')
        out_path = mydir + 'figs/SFS/D100/' + strain
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        for treatment in treatments:
            fig, ax = plt.subplots()
            reps = [x for x in IN.columns if 'frequency_L' + treatment in x]
            for rep in reps:
                sfs = IN[rep].values
                sfs = sfs[~np.isnan(sfs)]
                #weights = np.ones_like(sfs)/len(sfs)
                rep_number = rep[-1]
                if len(sfs) > 0:
                    weights = np.ones_like(sfs)/len(sfs)
                    ax.hist(sfs, 50, fc=colors[rep_number], histtype='stepfilled',
                        label='Replicate ' + rep_number, alpha=0.5, weights= weights)
            plt.xlim([0.01,0.99])
            ax.legend()
            ax.set_xlabel('Site frequency', fontsize = 16)
            ax.set_ylabel('Relative abundance', fontsize = 16)
            title = 'Site frequency spectra for ' + treatment_dict[rep[-3]] \
                + ' ' + species_dict[rep[-2]]
            fig.suptitle(title, fontsize=15)
            plt.savefig(out_path + '/' + strain + '_' + treatment + '.png', dpi=600)
            plt.close()

    #KDE = CV_KDE(afs)
    #fig, ax = plt.subplots()
    #ax.plot(KDE[0], KDE[1], linewidth=3, alpha=0.5, label='bw=%.2f' % KDE[2])
    #weights1 = np.ones_like(afs1)/len(afs1)
    #weights2 = np.ones_like(afs2)/len(afs2)
    #ax.hist(afs1, 50, fc='b', histtype='stepfilled', alpha=0.5)
    #ax.hist(afs2, 50, fc='r', histtype='stepfilled', alpha=0.5)
    #plt.xlim([0.01,0.99])
    #ax.set_xlabel('Site frequency', fontsize = 16  )
    #ax.set_ylabel('Number of sites', fontsize = 16)
    #title = 'Site frequency spectreum of '
    #ax.text(2, 1650, 'Site-specific substitutions', fontsize=15)
    #plt.savefig(mydir + 'figs/SFS/D100/D/test.png', dpi=600)
    #plt.close()

    #IN[np.isfinite(IN['frequency_L0D5'])]


#def SFS_figs():

def K_fig():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    K_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        K_values.append( IN[IN.strain == strain]['k_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], K_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)
    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_xticklabels(labels, fontsize = 18)
    ax.set_ylabel('Substitutions', fontsize = 20)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    #ax.set_xticklabels(treatments, fontsize = 8)
    #ax.get_xaxis().tick_bottom()
    #ax.get_yaxis().tick_left()
    #ax.set_ylim(1.15, 1.7)
    #ax.set_title('Genera')
    #ax.set_ylabel('Substitutions')
    #plt.plot([1, 1, 7, 7] , [1.6, 1.65, 1.65, 1.6], lw=1.5, c='k')
    #plt.text(4, 1.66, "****", ha='center', va='bottom', color='k', fontsize = 17)
    fig.savefig(mydir + 'figs/K.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def pi_fig():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    pi_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        pi_values.append( IN[IN.strain == strain]['S_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], pi_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)

    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_xticklabels(labels, fontsize = 18)
    ax.set_ylabel('Nucleotide diversity (pi)', fontsize = 20)
    plt.ticklabel_format(style='sci', axis='y')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    fig.savefig(mydir + 'figs/pi.png', bbox_inches='tight',  dpi = 600)
    plt.close()


def W_theta_fig():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    W_theta_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        W_theta_values.append( IN[IN.strain == strain]['W_T_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], W_theta_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)

    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_xticklabels(labels, fontsize = 18)
    ax.set_ylabel('Wattersons theta', fontsize = 20)
    plt.ticklabel_format(style='sci', axis='y')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    fig.savefig(mydir + 'figs/WT.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def poly_fig():
    IN = pd.read_csv(mydir + 'data/allele_freqs/polymorphism_number.txt', sep = '\t', header = 'infer')
    strains = list(IN.Strains.unique())
    treatments = list(IN.Treatments.unique())
    poly_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        poly_values.append( IN[IN.Strains == strain]['Polymorphism_number'].tolist())
        treat_values.append( IN[IN.Strains == strain]['Treatments'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], poly_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)
    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_xticklabels(labels, fontsize = 18)

    ax.set_ylabel('Number of polymorphisms', fontsize = 20)

    fig.savefig(mydir + 'figs/poly.png', bbox_inches='tight',  dpi = 600)
    plt.close()

#poly_fig()
#K_fig()
#pi_fig()
#W_theta_fig()

AFS()
