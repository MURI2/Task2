from __future__ import division
import os, re
import pandas as pd

mydir = os.path.expanduser("~/GitHub/Task2/Prod/")


IN = pd.read_csv(mydir + 'data/raw_data/test.csv', sep = ',')
print IN
