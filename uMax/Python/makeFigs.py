from __future__ import division
import os, re
import gomp as gp
import pandas as pd
import  matplotlib.pyplot as plt


mydir = os.path.expanduser("~/GitHub/Task2/uMax/")

def getTransferTime(x):
    if x[1] == str(0):
        return 1
    elif x[1] == str(1):
        return 10
    elif x[1] == str(2):
        return 100



path_IN = mydir + 'data/raw_data/Task2_24hr_24well_2017_05_25.txt'
path_OUT = mydir + 'data/clean_data/Task2_24hr_24well_2017_05_25.txt'

#gp.cleanData(path_IN, path_OUT, wells = 24)
#gp.modGompGrowth(path_OUT, smooth = True)

params_IN = mydir + 'data/params/Task2_24hr_24well_2017_05_25.txt'
IN = pd.read_csv(params_IN, sep = ' ')

IN['TransferTime'] = IN['Sample'].apply(getTransferTime)
umax_mean = IN['umax'].groupby(IN['Sample']).mean().values
A_mean = IN['A'].groupby(IN['Sample']).mean().values
L_mean = IN['L'].groupby(IN['Sample']).mean().values
x = IN['TransferTime'].groupby(IN['Sample']).mean().values
print IN['L'].groupby(IN['Sample']).mean()

fig = plt.figure()
plt.scatter(x, umax_mean, c='#87CEEB', marker='o', label='_nolegend_', s = 60)
plt.title('100 day Janthino', fontsize = 24)
plt.xlabel('Transfer time', fontsize = 18)
plt.ylabel('maximum growth rate', fontsize = 18)
plt.xscale('log')
fig_name = mydir + 'figs/umax.png'
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

fig = plt.figure()
plt.scatter(x, A_mean, c='#87CEEB', marker='o', label='_nolegend_', s = 60)
plt.title('100 day Janthino', fontsize = 24)
plt.xlabel('Transfer time', fontsize = 18)
plt.ylabel('Yield', fontsize = 18)
plt.xscale('log')
fig_name = mydir + 'figs/A.png'
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()


fig = plt.figure()
plt.scatter(x, L_mean, c='#87CEEB', marker='o', label='_nolegend_', s = 60)
plt.title('100 day Janthino', fontsize = 24)
plt.xlabel('Transfer time', fontsize = 18)
plt.ylabel('Lag parameter (lambda)', fontsize = 18)
plt.xscale('log')
fig_name = mydir + 'figs/L.png'
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
