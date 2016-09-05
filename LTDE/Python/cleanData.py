from __future__ import division
import pandas as pd
import os, re
import numpy as np

mydir = os.path.expanduser("~/github/Task2/LTDE")

strain_dict = {'ATCC13985': 'Pseudomonas', 'ATCC43928': 'Pseudomonas', \
'KBS0702': 'Arthrobacter', 'KBS0703': 'Arthrobacter', 'KBS0705':'Inquilinus', \
'KBS0706':'Mycobacterium', 'KBS0707': 'Pseudomonas', 'KBS0710': 'Pseudomonas', \
'KBS0711': 'Janthinobacterium', 'KBS0712': 'Variovorax', 'KBS0713': 'Yersinia', \
'KBS0715': 'Curtobacterium', 'KBS0721': 'Flavobacterium', 'KBS0722': 'Oerskovia', \
'KBS0724': 'Rhodococcus', 'KBS0727': 'Bradyrhizobium', 'KBS0801': 'Burkholderia', \
'KBS0802': 'Pseudomonas', 'KBS0812': 'Bacillus'}


def clean_iRep(path = mydir):
    #strains = ['KBS0703', 'KBS0710', 'KBS0711', 'KBS0713', \
    #    'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802', 'KBS0812']
    strains = ['ATCC13985', 'ATCC43928', 'KBS0702', 'KBS0703', 'KBS0705', \
        'KBS0706', 'KBS0707', 'KBS0710', 'KBS0711','KBS0712', 'KBS0713', \
        'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0801', \
        'KBS0802', 'KBS0812']
    lines_to_keep = [6, 10, 14, 18]
    to_pandas = []
    for strain in strains:
        strainPath = path + '/data/iRep/' + strain + '/'
        for content in os.listdir(strainPath):
            if content.endswith('.tsv'):
                # iRep, r2, coverage, % windows passing filer.
                name = re.split(r'[._]+', content)
                content_list = []
                content_list.append(name[0])
                content_list.append(name[1])
                for count, line in enumerate(open(strainPath + content)):
                    split_line = re.split(r'[#\t\n]+', line)
                    # remove empty items
                    clean_line = filter(None, split_line)
                    if len(clean_line) == 0:
                        continue
                    #print clean_line

                    if count in lines_to_keep:
                        if clean_line[-1] == 'n/a':
                            content_list.append(np.nan)
                        else:
                            content_list.append(float(clean_line[-1]))
                to_pandas.append( content_list)
    df = pd.DataFrame(to_pandas)
    df.columns  =['Strain', 'Replicate', 'iRep', 'R2', 'Coverage', 'PercentPass']
    df_path = path + '/data/iRep/iRepFinal.txt'
    df.to_csv(df_path, sep=' ', index = False)



def mergeData(path = mydir):
    num_to_letter = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', \
        '7': 'G', '8':'H', '9':'I', '10':'J', '11':'K', '12':'L', '13':'M'}
    IN_iRep = pd.read_csv(mydir + '/data/iRep/iRepFinal.txt', sep =' ' )
    IN_PopGen = pd.read_csv(mydir + '/data/mapgd/final/PopGenStats.txt', sep = ' ' )
    IN_merged_1 = pd.merge(IN_iRep, IN_PopGen, how='right', on=['Strain', 'Replicate'])
    pd.set_option('display.precision',20)
    IN_Evol = pd.read_csv(mydir + '/data/perRepDeathCurveTraits_centered.txt', sep = ' ')
    print IN_Evol
    IN_Evol['log_evolvability'] = np.log(IN_Evol.evolvability)
    replicate = []
    # model_type says whether the chosem model was linear (0) or quadratic (1)
    model_type = []
    for i, row in IN_Evol.iterrows():
        if row['evolvability'] == 0:
            model_type.append(0)
        else:
            model_type.append(1)

        if row['strain'] == 'KBS0711W':
            IN_Evol.set_value(i,'strain','KBS0711')
            if row['rep'] == 1:
                replicate.append('WA')
            elif row['rep'] == 2:
                replicate.append('WB')
            elif row['rep'] == 3:
                replicate.append('WC')
            elif row['rep'] == 4:
                replicate.append('WD')
        else:
            replicate.append(num_to_letter[str(row['rep'])])
    IN_Evol['Replicate'] = np.asarray(replicate)
    IN_Evol['model_type'] = np.asarray(model_type)
    IN_Evol.rename(columns={'strain':'Strain'}, inplace=True)
    IN_merged_2 = pd.merge(IN_merged_1, IN_Evol, how='left', on=['Strain', 'Replicate'])
    del IN_merged_2['rep']
    #IN_merged_2 = IN_merged_2[IN_merged_2.Strain != 'KBS0715']
    #IN_merged_2 = IN_merged_2[IN_merged_2.Strain != 'KBS0722']
    IN_merged_2 = IN_merged_2[IN_merged_2.T_D > 0]
    #IN_merged_2 = IN_merged_2[IN_merged_2.iRep <= 2.2]

    IN_merged_2.to_csv(path + '/data/Final/MAPGD_Evol_iRep.txt',sep='\t', index=False)

#clean_iRep()
mergeData()
