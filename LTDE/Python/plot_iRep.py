from __future__ import division
import pandas as pd
import os, math
import numpy as np
import csv, collections
from itertools import chain
from scipy import stats
import  matplotlib.pyplot as plt

mydir = os.path.expanduser("~/github/Task2/LTDE")

def plot_iRep():
    strains = ['KBS0703', 'KBS0705', 'KBS0706', 'KBS0710', 'KBS0711', 'KBS0713', \
        'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802']
    content_list = []
    for content in os.listdir(mydir + '/data/GATK/raw/' + strain): # "." means current directory
        content_list.append(content)
