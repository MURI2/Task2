from __future__ import division
import pandas as pd
import numpy as np
import os, re
from datetime import date
import  matplotlib.pyplot as plt

mydir = os.path.expanduser("~/github/Task2/LTDE/")

data = pd.read_csv(mydir + 'data/longtermdormancyAugust16.csv')
data = data.iloc[:,:-1]

#strains = ['KBS0722', 'KBS0724', 'KBS0727', 'KBS0715', 'KBS0703', 'KBS0711', \
#            'KBS0713', 'KBS0802']

strains = ['KBS0711']


def dateToDays(row):
    mont_dict = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6, 'Jul':7, \
        'Aug':8, 'Sep':9, 'Oct':10, 'Nov': 11, 'Dec':12}
    start = row['Dormstart_date']
    end = row['Firstread_date']
    start_split = re.split(r'[-_]+', start)
    end_split = re.split(r'[-_]+', end)

    #end.split('-')
    year_start = int('20' + start_split[2])
    year_end = int('20' + end_split[2])
    d0 = date(year_start, mont_dict[start_split[1]], int(start_split[0]))
    d1 = date(year_end, mont_dict[end_split[1]], int(end_split[0]))
    delta = d1 - d0
    return delta.days

def test_iter():
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for strain in strains:
        data_strain = data.loc[data['Strain'] == strain]
        data_strain_unique = data_strain.Rep.unique()
        for x in data_strain_unique:
            rep = data_strain.loc[data_strain['Rep'] == x]
            rep['Days'] = rep.apply(dateToDays, axis=1)
            days = rep['Days'].values
            colonies = rep['Colonies'].values
            zip_test = zip(days, colonies)
            zip_test_filter = [x for x in zip_test if x[1] != 'fungus']
            print zip_test
            ax1.scatter(zip_test_filter[0], zip_test_filter[1], color='blue', edgecolor='none')
    fig.savefig(mydir + '/figs/test_colonies.png',  bbox_inches = "tight", pad_inches = 0.4, dpi = 600)




#slope, intercept, r_value, p_value, std_err = stats.linregress(TD,slope)
#print slope, r_value, p_value
#plt.xlim([0,100])
#plt.xlabel('Tajimas D', fontsize=20)
#plt.ylabel('Slope', fontsize=20)

test_iter()
