from __future__ import division
import pandas as pd
import os, math
import numpy as np
import  matplotlib.pyplot as plt
import csv

mydir = os.path.expanduser("~/github/Task2/LTDE")

def plotTD():
    IN = (mydir + '/data/mapgd/final/mapgd_TD.txt')
    strains = []
    data = []
    with open(IN) as f:
        my_lines = f.readlines()
        for x in my_lines:
            x = x.strip().split(',')
            if x[0] == 'KBS0721' or x[0] == 'KBS0710':
                continue
            strains.append(x[0])
            data_x = [ float(y) for y in x[1:] ]
            print np.mean(data_x)
            data.append(data_x)

    fig = plt.figure(1, figsize=(9, 6))

    # Create an axes instance
    ax = fig.add_subplot(111)
    ax.set_ylim(-3.5, 3.5)
    ax.axhline(linewidth=2, color='darkgrey',ls='--')


    # Create the boxplot
    bp = ax.boxplot(data)
    ax.set_xticklabels(strains)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_ylabel('Tajimas D')
    # Save the figure
    fig.savefig(mydir + '/figs/Tajimas_D.png', bbox_inches='tight')

plotTD()
