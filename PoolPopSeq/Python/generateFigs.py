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
import sklearn.metrics.pairwise as pairwise
import skbio.stats.ordination as ordination
import skbio.stats.distance as distance
from matplotlib.patches import Polygon

import pylab as P


mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/")

strain_colors = {'Janthinobacterium':'indigo', 'Caulobacter':'lightblue', \
         'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'cyan',
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
                sfs = sfs[(sfs != float(1))]
                #weights = np.ones_like(sfs)/len(sfs)
                rep_number = rep[-1]
                if len(sfs) > 0:
                    weights = np.ones_like(sfs)/len(sfs)
                    ax.hist(sfs, 50, fc=colors[rep_number], histtype='stepfilled',
                        label='Replicate ' + rep_number, alpha=0.5, weights= weights)
            plt.xlim([0.01,0.99])
            ax.legend()
            ax.set_xlabel('Site frequency', fontsize = 16)
            ax.set_ylabel('Fraction of sites', fontsize = 16)
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
    #K_values = np.log10(K_values)
    for i, color in enumerate(colors):
        ax.plot(treat_values[i], K_values[i], marker='o', linestyle='', \
            ms=6, label=strains[i], color = color, alpha = 0.9)
    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[1] = '1-Day'
    labels[3] = '10-Day'
    labels[5] = '100-Day'
    ax.set_yscale("log", nonposy='clip')
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

def poly_fig(variable):
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    poly_values = []
    treat_values = []
    strain_values = []
    colors = []
    strain_values_colors = []
    for treatment in treatments:
        for strain in strains:
            poly_values.append( IN[(IN.strain == strain) & (IN.treatment == treatment)][variable].tolist() )
            #treat_values.append( IN[(IN.strain == strain) & (IN.treatment == treatment)]['treatment'].tolist() )
            #strain_values.append( IN[(IN.strain == strain) & (IN.treatment == treatment)]['strain'].tolist() )
            treat_values.append(treatment)
            strain_values.append(strain)
            strain_values_colors.append(strain_colors[strain])
        colors.append(strain_colors[strain])
    poly_values_per_gen = []
    for i, poly_value in enumerate(poly_values):
        if variable == 'W_T_L':
            poly_value = [x for x in poly_value if x != 0]
        if treat_values[i] == 0:
            poly_values_per_gen.append([ (x / 100) for x in poly_value])
        elif treat_values[i] == 1:
            poly_values_per_gen.append([ (x / 10) for x in poly_value])
        elif treat_values[i] == 2:
            poly_values_per_gen.append(poly_value)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('A Boxplot Example')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = plt.boxplot(poly_values_per_gen, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')

    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_xlabel('Transfer time', fontsize = 20)
    if variable == 'W_T_L':
        ax1.set_ylabel('Segregating sites, ' +  r'$\theta_{W}$' +  \
        ' \n per base per transfer', fontsize = 20)
    elif variable == 'pi_L':
        ax1.set_ylabel('Nucleotide diversity, ' +  r'$\pi$' +  \
        '\n per base per transfer', fontsize = 20)

    # Now fill the boxes with desired colors
    # numBoxes = number treatments * number taxa
    numBoxes = 3*6
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=strain_values_colors[i])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot([np.average(med.get_xdata())], [np.average(poly_values_per_gen[i])],
                 color='w', marker='*', markeredgecolor='k')
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels_new = []
    for i, label in enumerate(labels):
        if i == 3:
            labels_new.append('1-Day')
        elif i == 8:
            labels_new.append('10-Day')
        elif i == 14:
            labels_new.append('100-Day')
        else:
            labels_new.append('')
    ax1.set_xticklabels(labels_new, fontsize = 18)

    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)
    #strain_colors = {'Janthinobacterium':'cyan', 'Caulobacter':'lightblue', \
    #     'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'indigo',
    #     'Pedobacter': 'darkred'}
    plt.figtext(0.10, 0.640, 'Janthinobacterium',
                backgroundcolor='indigo',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.685,' Caulobacter',
            backgroundcolor='lightblue', color='black', weight='roman',
            size='x-small')
    plt.figtext(0.10, 0.730, 'Bacillus',
                backgroundcolor='cyan',
                color='black', weight='roman', size='x-small')
    plt.figtext(0.10, 0.775, 'Deinococcus',
                backgroundcolor='red',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.820, 'Pseudomonas',
                backgroundcolor='darkgreen',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.865, 'Pedobacter',
                backgroundcolor='darkred',
                color='white', weight='roman', size='x-small')

    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)


    ax1.set_yscale("log", nonposy='clip')
    fig.savefig(mydir + 'figs/' + variable +'.png', bbox_inches='tight',  dpi = 600)
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
        pi_values.append( IN[IN.strain == strain]['pi_L'].tolist())
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

def T_D_fig():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', \
        sep = ' ', header = 'infer')
    #IN = IN.sort(['treatment', 'strain'], ascending=[True, True])
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    T_D_values = []
    treat_values = []
    strain_values = []
    colors = []
    strain_values_colors = []
    for treatment in treatments:
        for strain in strains:
            T_D_value =  IN[(IN.strain == strain) & (IN.treatment == treatment)]['T_D'].tolist()
            T_D_values.append([x for x in T_D_value if np.isnan(x) == False])
            treat_values.append(treatment)
            strain_values.append(strain)
            strain_values_colors.append(strain_colors[strain])
        colors.append(strain_colors[strain])


    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('A Boxplot Example')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    bp = plt.boxplot(T_D_values, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_xlabel('Transfer time', fontsize = 20)
    ax1.set_ylabel(r'$D_{T }$', fontsize = 26)

    # Now fill the boxes with desired colors
    # numBoxes = number treatments * number taxa
    numBoxes = 3*6
    medians = list(range(numBoxes))
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=strain_values_colors[i])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot([np.average(med.get_xdata())], [np.average(T_D_values[i])],
                 color='w', marker='*', markeredgecolor='k')
    #labels = [item.get_text() for item in ax1.get_xticklabels()]
    #labels_new = []
    #for i, label in enumerate(labels):
    #    if i % 3 == 0:
    #        labels_new.append('1-Day')
    #    elif i % 3 == 1:
    #        labels_new.append('10-Day')
    #    elif i % 3 == 2:
    #        labels_new.append('100-Day')

    #ax1.set_xticklabels(labels_new, fontsize = 18)

    #for label in ax1.get_xmajorticklabels():
    #    label.set_rotation(60)
    #    label.set_fontsize(8)
        #print ', '.join(i for i in dir(label) if not i.startswith('__'))
    #    label.set_horizontalalignment("right")

    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels_new = []
    for i, label in enumerate(labels):
        if i == 3:
            labels_new.append('1-Day')
        elif i == 8:
            labels_new.append('10-Day')
        elif i == 14:
            labels_new.append('100-Day')
        else:
            labels_new.append('')
    ax1.set_xticklabels(labels_new, fontsize = 18)

    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)
    #strain_colors = {'Janthinobacterium':'cyan', 'Caulobacter':'lightblue', \
    #     'Deinococcus': 'red', 'Pseudomonas':'darkgreen', 'Bacillus':'indigo',
    #     'Pedobacter': 'darkred'}
    plt.figtext(0.10, 0.640, 'Janthinobacterium',
                backgroundcolor='indigo',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.685,' Caulobacter',
            backgroundcolor='lightblue', color='black', weight='roman',
            size='x-small')
    plt.figtext(0.10, 0.730, 'Bacillus',
                backgroundcolor='cyan',
                color='black', weight='roman', size='x-small')
    plt.figtext(0.10, 0.775, 'Deinococcus',
                backgroundcolor='red',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.820, 'Pseudomonas',
                backgroundcolor='darkgreen',
                color='white', weight='roman', size='x-small')
    plt.figtext(0.10, 0.865, 'Pedobacter',
                backgroundcolor='darkred',
                color='white', weight='roman', size='x-small')


    fig.savefig(mydir + 'figs/T_D.png', bbox_inches='tight',  dpi = 600)
    plt.close()





def T_D_figsss():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', \
        sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    W_theta_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        print IN[IN.strain == strain]['T_D'].tolist()
        W_theta_values.append( IN[IN.strain == strain]['T_D'].tolist())
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
    ax.set_ylabel('Tajimas D', fontsize = 20)
    #plt.ticklabel_format(style='sci', axis='y')
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    fig.savefig(mydir + 'figs/TD.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def pi_vs_k2():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    pi_values = []
    k_values = []
    treat_values = []
    strain_values = []
    colors = []
    for strain in strains:
        pi_values.append( IN[IN.strain == strain]['pi_L'].tolist())
        k_values.append( IN[IN.strain == strain]['k_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for i, color in enumerate(colors):
        #ax.scatter(k_values[i], pi_values[i], marker='o', linestyle='', \
        #    ms=6, label=strains[i], color = color, alpha = 0.9)
        print k_values[i]
        print pi_values[i]
        ax.scatter(k_values[i], pi_values[i], marker='o', \
            label=strains[i], color = color, alpha = 0.9)

    ax.legend(numpoints=1, prop={'size':14},  loc='upper right', frameon=False)
    fig.canvas.draw()

    #labels = [item.get_text() for item in ax.get_xticklabels()]
    #labels[1] = '1-Day'
    #labels[3] = '10-Day'
    #labels[5] = '100-Day'
    ax.set_ylabel('Substitutions per-site (K)', fontsize = 20)
    ax.set_ylabel('Nucleotide diversity per-site (pi)', fontsize = 20)
    #ax.set_xscale("log", nonposy='clip')
    #ax.set_yscale("log", nonposy='clip')
    plt.ticklabel_format(style='sci', axis='y')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    fig.savefig(mydir + 'figs/pi_v_k.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def pi_vs_k():
    IN = pd.read_csv(mydir + 'data/pop_gen_stats/D100/popGenTable.txt', sep = ' ', header = 'infer')
    strains = list(IN.strain.unique())
    treatments = list(IN.treatment.unique())
    pi_values = []
    k_values = []
    treat_values = []
    strain_values = []
    colors = []
    print IN[IN.strain == 'Pedobacter']
    for strain in strains:
        pi_values.append( IN[IN.strain == strain]['pi_L'].tolist())
        k_values.append( IN[IN.strain == strain]['k_L'].tolist())
        treat_values.append( IN[IN.strain == strain]['treatment'].tolist())
        colors.append(strain_colors[strain])
        #strain_values.append(IN[IN.Strains == strain]['Strains'].tolist())
    fig, ax = plt.subplots()
    for i, color in enumerate(colors):
        ax.plot(k_values[i], pi_values[i], marker='o', alpha = 0.8, \
            linestyle='', ms=12, c = color)
        #ax.scatter(k_values[i], pi_values[i], marker='o', \
        #    label=strains[i], color = color, alpha = 0.9)
    ax.legend(numpoints=1, prop={'size':10},  loc='upper left', frameon=False)
    #slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    #print "slope = " + str(slope)
    #print "r2 = " + str(r_value**2)
    #print "p = " + str(p_value)
    #predict_y = intercept + slope * x
    #pred_error = y - predict_y
    #degrees_of_freedom = len(x) - 2
    #residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    ax.set_xlabel('Substitutions per-site (K)', fontsize = 20)
    ax.set_ylabel('Nucleotide diversity per-site (pi)', fontsize = 20)
    ax.set_xscale("log", nonposy='clip')
    sax.set_yscale("log", nonposy='clip')

    fig.savefig(mydir + 'figs/pi_v_k.png', bbox_inches='tight',  dpi = 600)
    plt.close()

def sample_by_gene_dissimilarity(day, strain):
    path = mydir + 'data/gene_by_sample/' + strain + '/' + \
        day + '/sample_by_gene_multiple_hits.txt'
    IN = pd.read_csv(path, sep = '\t', header = 'infer', index_col = 0)
    if day == 'D100' and strain == 'B':
        IN = IN.drop(['frequency_L2B3', 'frequency_L1B1'], axis=0)
    IN = IN[(IN.T != 0).any()]
    rows = IN.index
    IN_np = IN.as_matrix()
    distance_sklearn = pairwise.pairwise_distances(IN_np, metric='braycurtis')
    distance_skbio = distance.DistanceMatrix(distance_sklearn, ids=rows)
    pcoa = ordination.PCoA(distance_skbio)
    pcoa_results = pcoa.scores()
    #print pcoa_results.__dict__
    names = [i.split('_')[1] for i in pcoa_results.site_ids]
    x = pcoa_results.site[:,0]
    y = pcoa_results.site[:,1]
    zipped = zip(x, y, names)
    fig, ax = plt.subplots()
    treatments = [0, 1, 2]
    treatment_names = ['1-Day', '10-Day', '100-Day']
    treatment_colors = ['#87CEEB', '#FFA500', '#FF6347']
    for treatment in treatments:
        zipped_treatment = [k for k in zipped if str(treatment) in k[2]]
        x_treatment = [k[0] for k in zipped_treatment]
        y_treatment = [k[1] for k in zipped_treatment]
        ax.plot(x_treatment, y_treatment, marker='o', alpha = 0.8, \
            linestyle='', ms=12, label=treatment_names[treatment], \
            c = treatment_colors[treatment])
    for item in zipped:
        ax.annotate(item[2], (item[0],item[1]))
        round(14.22222223, 2)
    xlabel = 'PCoA 1 (' +  str(round(pcoa_results.proportion_explained[0] * 100, 1))  + '%)'
    ylabel = 'PCoA 2 (' +  str(round(pcoa_results.proportion_explained[1] * 100, 1))  + '%)'
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    title = 'Day 100 ' + species_dict[strain]
    plt.title(title, fontsize = 22)
    fig.savefig(mydir + 'figs/pcoa/D100_' + strain + '.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def GMD_hist(day, strain):
    path = mydir + 'data/gene_by_sample/' + strain + '/' + day + \
        '/sample_by_gene_multiple_hits.txt'
    IN = pd.read_csv(path, sep = '\t', header = 'infer', index_col = 0)
    sample_names = IN.index.values
    day_1 = [x for x in sample_names if x[-3] == str(0) ]
    day_10 = [x for x in sample_names if x[-3] == str(1) ]
    day_100 = [x for x in sample_names if x[-3] == str(2) ]
    day_1_sum = IN.loc[day_1,:].sum(axis=0)
    day_10_sum = IN.loc[day_10,:].sum(axis=0)
    day_100_sum = IN.loc[day_100,:].sum(axis=0)
    day_1_sum.name = '1-Day'
    day_10_sum.name = '10-Day'
    day_100_sum.name = '100-Day'

    merged = pd.concat([day_1_sum, day_10_sum, day_100_sum], axis=1)
    idx = merged.sum(axis=1).sort_values(ascending=False).head(30).index
    merged_sort_descending = merged.ix[idx]
    gene_names = merged_sort_descending.index.values

    fig, ax = plt.subplots()
    bar_width = 0.65
    # positions of the left bar-boundaries
    bar_l = [i+1 for i in range(len(merged_sort_descending['100-Day']))]
    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i+(bar_width/2) for i in bar_l]
    # Create a bar plot, in position bar_1
    ax.bar(bar_l,
            # using the 100-day data
            merged_sort_descending['100-Day'],
            # set the width
            width=bar_width,
            # with the label pre score
            label='100-Day',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#FF6347')

    # Create a bar plot, in position bar_1
    ax.bar(bar_l,
            # using the 10-day data
            merged_sort_descending['10-Day'],
            # set the width
            width=bar_width,
            # with 100-day on the bottom
            bottom=merged_sort_descending['100-Day'],
            # with the label mid score
            label='10-Day',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#FFA500')

    # Create a bar plot, in position bar_1
    ax.bar(bar_l,
            # using the 1-day data
            merged_sort_descending['1-Day'],
            # set the width
            width=bar_width,
            # with pre_score and mid_score on the bottom
            bottom=[i+j for i,j in zip(merged_sort_descending['100-Day'],\
                merged_sort_descending['10-Day'])],
            # with the label post score
            label='1-Day',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#87CEEB')

    # set the x ticks with names
    plt.xticks(tick_pos, gene_names)
    ax.set_ylabel("Number of mutations", fontsize = 16 )
    #ax.set_xlabel("Gene names")
    plt.legend(loc='upper right')
    for label in ax.get_xmajorticklabels():
        label.set_rotation(60)
        label.set_fontsize(8)
        #print ', '.join(i for i in dir(label) if not i.startswith('__'))
        label.set_horizontalalignment("right")
    # make room for the labels
    plt.gcf().subplots_adjust(bottom=0.20)

    # Set a buffer around the edge
    plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])

    title = 'Day ' + day[1:] + ' ' + species_dict[strain]
    plt.title(title, fontsize = 22)
    fig_dir = mydir + 'figs/GMD/' + day
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_out = fig_dir + '/' + strain + '.png'
    #plt.savefig(fig_out, dpi=600)
    fig.savefig(fig_out, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

#poly_fig()
#K_fig()
#pi_fig()
#W_theta_fig()
#T_D_fig()
#pi_vs_k()
#AFS()
#poly_fig('pi_L')
#poly_fig('W_T_L')
poly_fig('k_L')

#sample_by_gene_dissimilarity('D100', 'D')
#sample_by_gene_dissimilarity('D100', 'F')
#sample_by_gene_dissimilarity('D100', 'P')
#sample_by_gene_dissimilarity('D100', 'C')
#sample_by_gene_dissimilarity('D100', 'B')
#sample_by_gene_dissimilarity('D100', 'J')

#GMD_hist('D100', 'D')
#GMD_hist('D100', 'C')
#GMD_hist('D100', 'F')
#GMD_hist('D100', 'P')
#GMD_hist('D100', 'B')
#GMD_hist('D100', 'J')
