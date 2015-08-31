import numpy as np
import itertools
import random

shakers = ['T', 'B']
strains = ['A', 'B', 'C', 'D', 'F', 'J']
transfers = ['0', '1', '2']
replicates = 5

total = len(shakers) * len(strains) * len(transfers) * replicates

wells = np.asarray(['%s%02d' % (r, c) for r in 'ABCDEFGH' for c in range(1, 13)])

# Now split 96 wells into the blocks we want
#print wells[0:6] + wells[12:18]
B1 = np.concatenate([wells[0:6], wells[12:18], wells[24:30]])
B2 = np.concatenate([wells[6:12], wells[18:24], wells[30:36]])
B3 = np.concatenate([wells[60:66], wells[72:78], wells[84:90]])
B4 = np.concatenate([wells[66:72], wells[78:84], wells[90:96]])
B5 = np.concatenate([wells[38:47], wells[50:59]])
Blocks = [B1, B2, B3, B4, B5]

# Get all permutations of strains and transfer times
PermTrt = np.asarray([''.join(i) for i in list(itertools.product(strains, transfers))])

#print random.sample(zip(PermTrt,B1), len(B1))

for x in shakers:
    #structWells = np.reshape(wells, (-1, 12))
    structWells = wells.copy()
    count = 1
    for y in Blocks:
        #trtLaidDown = random.sample(zip(PermTrt,y), len(y))
        trtLaidDown1 = random.sample(PermTrt, len(PermTrt))
        trtLaidDown2 = random.sample(y, len(y))
        trtLaidDown = zip(trtLaidDown1, trtLaidDown2)
        for z in trtLaidDown:
            splitFirst = list(z[0])
            label = str(splitFirst[1]) + "-" + str(x) + "-" + str(z[1]) + "-" + str(splitFirst[0]) + "-" + str(count)
            np.where((structWells==str(z[1])),label,structWells)
            print label
        count += 1
