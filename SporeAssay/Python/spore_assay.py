from __future__ import division
import os
import pandas as pd
import  matplotlib.pyplot as plt

mydir = os.path.expanduser("~/GitHub/Task2/SporeAssay/")

def make_fig():
    path_IN = mydir + 'data/Sporulation_170822_long.txt'
    IN = pd.read_csv(path_IN, sep = '\t')
    IN_0B1 = IN.loc[IN['Pop'] == '0B1']
    IN_2B1 = IN.loc[IN['Pop'] == '2B1']
    IN_mean_0B1 = IN_0B1['Vegetative_percent'].groupby(IN_0B1['Time_hours']).mean().reset_index()
    IN_mean_2B1 = IN_2B1['Vegetative_percent'].groupby(IN_2B1['Time_hours']).mean().reset_index()

    IN_std_0B1 = IN_0B1['Vegetative_percent'].groupby(IN_0B1['Time_hours']).std().reset_index()
    IN_std_2B1 = IN_2B1['Vegetative_percent'].groupby(IN_2B1['Time_hours']).std().reset_index()

    fig = plt.figure()

    #plt.scatter(IN_mean_0B1.Time_hours.values, IN_mean_0B1.Vegetative_percent.values, c='#87CEEB', marker='o', label='_nolegend_', s = 60)
    plt.plot(IN_mean_0B1.Time_hours.values, 1.001- IN_mean_0B1.Vegetative_percent.values, \
        'b-',  c='#87CEEB')
    plt.plot(IN_mean_2B1.Time_hours.values, 1.001- IN_mean_2B1.Vegetative_percent.values, \
        'b-',  c = '#FF6347')
    plt.errorbar(IN_mean_0B1.Time_hours.values, 1.001- IN_mean_0B1.Vegetative_percent.values, \
        IN_std_0B1.Vegetative_percent.values,  linestyle='None', marker='o', c='#87CEEB', elinewidth=1.5, label="Bacillus population transferred daily",)
    plt.errorbar(IN_mean_2B1.Time_hours.values, 1.001- IN_mean_2B1.Vegetative_percent.values, \
        IN_std_2B1.Vegetative_percent.values, linestyle='None', marker='o', c = '#FF6347', elinewidth=1.5, label="Bacillus population transferred every 100 days",)
    #plt.title('Bacillus sporulation', fontsize = 24)
    plt.xlabel('Time (hrs)', fontsize = 18)
    plt.ylabel('Percent spores, log10', fontsize = 18)
    plt.ylim(0.0008, 1.1)
    plt.yscale('log')
    plt.legend(numpoints=1, prop={'size':10},  loc='lower right', frameon=False)
    fig_name = mydir + 'figs/test.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





make_fig()
