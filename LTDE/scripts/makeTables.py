from __future__ import division
import pandas as pd
import os

mydir = os.path.expanduser("~/github/Task2/LTDE")

# ignore 705 and 706

species = ['KBS0703', 'KBS0710', 'KBS0711', 'KBS0713', 'KBS0715', 'KBS0721', \
            'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802']


def pandasNames(dataFrame):
    names = ['Scaffold', 'Pos', 'Major', 'Minor', 'Coverage', 'Error', \
        'Gene', 'AA-change']
    samples = len(df.columns) - len(names)
    for x in range(samples):
        name = 'Sample_' + str(x+1)
        names.append(name)
    dataFrame.columns = names
    return dataFrame


test = mydir + '/data/mapgd/annotate/KBS0722_merged_annotate.pol'
df = pd.read_csv(test,  sep='\t', header = None)

df = pandasNames(df)
#df = df.drop('Sample_1', 1)
dfCoding = df.loc[df['AA-change'] != 'NC']
#print dfCoding
test = dfCoding.to_latex()
f = open('test.tex', 'w')
f.write(test.encode('utf8'))
f.close()
