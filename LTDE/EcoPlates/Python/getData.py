from __future__ import division
import pandas as pd
import os
# read it in, reshape to a 3 x 32 matrix and spit out as a .txt file
mydir = os.path.expanduser("~/github/Task2/LTDE/EcoPlates/")

def rawToLabelles(strains):
    resources  = ["Water", "Pyruvic.Acid.Methyl.Ester", "Tween.40", \
                       "Tween.80", "alpha-Cyclodextrin", "Glycogen", \
                       "D-Cellulobiose", "alpha-D-Lactose", "beta-Methyl-D-Glucoside", \
                       "D-Xylose", "i-Erythritol", "D-Mannitol", "N-Acetyl-D-Glucosamine", \
                       "D-Glucosaminic.Acid", "Glucose-1-Phosphate", \
                       "D,L-alpha-Glycerol.Phosphate", "D-Galactonic.Acid.gamma-Lactone", \
                       "D-Galacturonic.Acid", "2-Hydroxy.Benzoic.Acid", \
                       "4-Hydroxy.Benzoic.Acid", "gamma-Hydroxybutyric.Acid", \
                       "Itaconic.Acid", "alpha-Ketobutyric.Acid", "D-Malic.Acid", \
                       "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", \
                       "L-Threonine", "Glycyl-L-Glutamic.Acid", "Phenylethylamine", \
                       "Putrescine"]
    for strain in strains:
        IN_path = mydir + 'data/raw/' + strain
        OUT_path =  mydir + 'data/labelled/' + strain
        if not os.path.exists(OUT_path):
            os.makedirs(OUT_path)
        for xlsx in os.listdir(IN_path):
            if xlsx.endswith(".xlsx"):
                name = xlsx.split('.')[0]
                IN = IN_path + '/' + xlsx
                IN_read = pd.read_excel(IN, sheetname=0)
                shape = IN_read.shape
                if shape[0] == 8:
                    IN_data = IN_read.iloc[:,:-1]
                else:
                    IN_data = IN_read.iloc[4:,2:-1]
                print IN_data.shape
                arrays = []
                for x in range (0,12,4):
                    arrays.append(IN_data.iloc[:,x:x+4].values.flatten())
                data_labelled =  pd.DataFrame(arrays, columns = resources)
                OUT = OUT_path + '/' + name + '.txt'
                data_labelled.to_csv(OUT, sep='\t', index = False)




def normalizeData(strains):
    'Removes background using water and normalizes the data'
    for strain in strains:
        IN_path = mydir + 'data/labelled/' + strain
        OUT_path =  mydir + 'data/cleaned/' + strain
        if not os.path.exists(OUT_path):
            os.makedirs(OUT_path)
        for txt in os.listdir(IN_path):
            if txt.endswith(".txt"):
                name = txt.split('.')[0]
                name = name + '_cleaned'
                IN = IN_path + '/' + txt
                IN_read = pd.read_csv(IN, sep = '\t', index_col= False)
                IN_read_sub = IN_read.sub(IN_read.iloc[:,0].mean(),axis=0)
                # remove water column
                IN_read_sub_noW = IN_read_sub.iloc[:,1:]
                # negative values to zero
                IN_read_sub_noW[IN_read_sub_noW < 0] = 0
                IN_read_sub_noW = IN_read_sub_noW.round(4)
                OUT = OUT_path + '/' + name + '.txt'
                IN_read_sub_noW.to_csv(OUT, sep='\t', index = False)

strains = ['KBS0711']
rawToLabelles(strains)
normalizeData(strains)
