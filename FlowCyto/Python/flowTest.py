from __future__ import division
import FlowCytometryTools
from FlowCytometryTools import test_data_dir, test_data_file, FCMeasurement, FCPlate, ThresholdGate, PolyGate
import os, re, datetime, signal, collections
import matplotlib.pyplot as plt
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
from operator import itemgetter
import numpy as np
import pandas as pd
import scipy.stats as stats

data_path = os.path.expanduser('~/Box Sync/Flow_Cyto_Task2/Day100/')
git_path = os.path.expanduser('~/github/Task2/FlowCyto/')

# FITC-A = RSG
# Pacific Blue-A = DAPI
# APC-A = eFluor

class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException

def CV_KDE(oneD_array):
    # remove +/- inf
    oneD_array = oneD_array[np.logical_not(np.isnan(oneD_array))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(oneD_array[:, None])
    x_grid = np.linspace(np.amin(oneD_array), np.amax(oneD_array), 10000)
    kde = grid.best_estimator_
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple

def prepFileNames():
    strains = ['B', 'C', 'D', 'F', 'P', 'J']
    transfers = ['D1', 'D10', 'D100']
    reps = ['1', '2', '3', '4', '5']
    rename_dict = {'A':'A1', 'B':'A2', 'C':'A3', 'D':'A4', 'E':'A5', 'F':'A6',\
        'G':'A7', 'H':'A8', 'I':'A9', 'J':'A10'}
    for strain in strains:
        for transfer in transfers:
            for rep in reps:
                rep_path = data_path + strain + '/' + transfer + '/' + rep + '/'
                for filename in os.listdir(rep_path):
                    if filename.startswith('Well'):
                        pass
                    else:
                        filename_split = re.split(r'[_.]+', filename)
                        if filename_split[-2] in rename_dict.keys():
                            well_name = 'Well_' +  rename_dict[filename_split[-2] ]+ '_' + filename
                            os.rename(rep_path + filename, rep_path + well_name)


def makeFigFolders(analysis = 'RSG_DAPI'):
    strains = ['B', 'C', 'D', 'F', 'P', 'J']
    transfers = ['D1', 'D10', 'D100']
    reps = ['1', '2', '3', '4', '5']
    if analysis == 'RSG_DAPI':
        for strain in strains:
            for transfer in transfers:
                for rep in reps:
                    rep_path = git_path + 'figs/' + analysis + '/' + strain + '/' + transfer + '/' + rep + '/'
                    if not os.path.exists(rep_path):
                        os.makedirs(rep_path)
    elif analysis == 'RSG':
        for strain in strains:
            for transfer in transfers:
                path = git_path + 'figs/' + analysis + '/' + strain + '/' + transfer + '/'
                if not os.path.exists(path):
                    os.makedirs(path)
    #else:
    #    continue


def getDAPIgate(plate):
    As = ['A3', 'A4']
    cutoffs = []
    for A in As:
        DAPI = plate[A].data[['Pacific Blue-A']].values
        DAPI_gate = np.mean(DAPI) + (2*np.std(DAPI))
        cutoffs.append(DAPI_gate)
    cutoff = np.mean(cutoffs)
    return cutoff

def get_DAPI_RGF_Figs():
    strains = ['B', 'C', 'D', 'F', 'P', 'J']
    transfers = ['D1', 'D10', 'D100']
    reps = ['1', '2', '3', '4', '5']
    As = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8']
    for strain in strains:
        for transfer in transfers:
            for rep in reps:
                print strain, transfer, rep
                path = data_path + strain + '/' + transfer + '/' + rep + '/'
                plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
                plate = plate.dropna()
                plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
                    'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
                threshold = getDAPIgate(plate)
                for A in As:
                    fig = plt.figure()
                    DAPI = plate[A].data[['Pacific Blue-A']].values
                    RSG = plate[A].data[['FITC-A']].values
                    plt.scatter(DAPI, RSG)
                    plt.axvline(x = threshold)
                    name =  'RSG_DAPI' + strain + '_' + transfer + '_' + rep + '_' + A
                    plt.savefig(git_path + 'figs/RSG_DAPI/' + strain + '/' + transfer + \
                        '/' + rep + '/' + name + '.png')

def getDistRSG(KDEs = False):
    makeFigFolders(analysis = 'RSG')
    strains = ['B', 'C', 'D', 'F', 'P', 'J']
    transfers = ['D1', 'D10', 'D100']
    reps = ['1', '2', '3', '4', '5']
    #strains = ['C']
    #transfers = ['D100']
    #reps = ['4']
    OUT1 = open(git_path +'data/bandwidth.txt', 'w+')
    if KDEs == True:
        print>> OUT1, 'strain', 'transfer', 'rep', 'bandwidth', 'mean', 'std', 'median', 'skew'
        for strain in strains:
            for transfer in transfers:
                for rep in reps:
                    signal.signal(signal.SIGALRM, timeout_handler)
                    print strain, transfer, rep
                    path = data_path + strain + '/' + transfer + '/' + rep + '/'
                    plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
                    plate = plate.dropna()
                    plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
                        'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
                    threshold = getDAPIgate(plate)
                    gate = ThresholdGate(threshold, 'Pacific Blue-A', region='above')
                    gated_sample = plate['A8'].gate(gate)
                    RSG = gated_sample.data[['FITC-A']].values
                    signal.alarm(500)
                    a = datetime.datetime.now()
                    try:
                        KDE = CV_KDE(RSG)

                        df = pd.DataFrame({'x_grid':KDE[0],'pdf':KDE[1]})
                        if not os.path.exists(git_path + 'data/RSG'):
                            os.makedirs(git_path + 'data/RSG')
                        bw_round = round(KDE[2], 3)
                        pandas_name = 'RSG_'+ strain + '_' + transfer + '_' + rep + '.txt'
                        df.to_csv(path_or_buf= git_path + 'data/RSG/' + pandas_name, sep = ' ', index=False)
                        fig, ax = plt.subplots()
                        ax.plot(KDE[0], KDE[1], linewidth=3, alpha=0.5, label='bw=%.2f' % KDE[2])
                        ax.hist(RSG, 30, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
                        ax.legend(loc='upper left')
                        name = 'RSG_' + strain + '_' + transfer + '_' + rep
                        plt.savefig(git_path + 'figs/RSG'  + '/' + strain + '/' + transfer + '/' + name + '.png')
                        plt.close()

                        print>> OUT1, strain, transfer[1:], rep, str(KDE[2]), str(np.mean(RSG)), str(np.std(RSG)), \
                            str(np.median(RSG)), str(stats.skew(RSG)[0])
                        b = datetime.datetime.now()
                        c = b - a
                        print str(c.seconds) + " seconds"

                    except TimeoutException:
                        continue # continue the for loop if function takes more than x seconds
                    else:
                        # Reset the alarm
                        signal.alarm(0)
    else:
        print>> OUT1, 'Strain', 'Transfer', 'Rep', 'Mean', 'Std', 'Median', 'Skew'
        for strain in strains:
            for transfer in transfers:
                for rep in reps:
                    print strain, transfer, rep
                    path = data_path + strain + '/' + transfer + '/' + rep + '/'
                    plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
                    plate = plate.dropna()
                    plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
                        'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
                    threshold = getDAPIgate(plate)
                    gate = ThresholdGate(threshold, 'Pacific Blue-A', region='above')
                    gated_sample = plate['A8'].gate(gate)
                    RSG = gated_sample.data[['FITC-A']].values
                    print>> OUT1, strain, transfer[1:], rep, str(np.mean(RSG)), str(np.std(RSG)), \
                        str(np.median(RSG)), str(stats.skew(RSG)[0])



    OUT1.close()


def plotRSG():
    strain_dict = {'A':['cyan', 'Arthrobacter'], 'B':['lightblue', 'Bacillus'], \
            'C':['lightgreen', 'Caulobacter'], 'D':['darkred', 'Deinococcus'], 'F':['orange', 'Pedobacter'], \
            'J':['indigo', 'Janthinobacterium'], 'P':['darkgreen', 'Pseudomonas']}
    IN = pd.read_csv(git_path +'data/bandwidth.txt', sep = ' ')
    IN = IN.sort_values(by=["Transfer","Skew"], ascending=[True, False])
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))
    x = np.log10(IN.Transfer.values)
    y = IN.Skew.values
    groups = IN.groupby('Strain')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(np.log10(group.Transfer), group.Skew, marker='o', alpha = 0.4 , \
            linestyle='', ms=12, label=strain_dict[name][1], c = strain_dict[name][0])
    ax.legend(numpoints=1, prop={'size':10},  loc='upper left', frameon=False)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print "slope = " + str(slope)
    print "r2 = " + str(r_value**2)
    print "p = " + str(p_value)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    degrees_of_freedom = len(x) - 2
    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
    plt.plot(x, predict_y, 'k-')
    plt.tick_params(axis='both', which='major', labelsize=10)
    #plt.ylim([min(y)-0.1,max(y)+0.1])
    #plt.xlim([min(x)-0.001,max(x)+0.001])
    plt.xlabel('Days starved, log10', fontsize=20)
    plt.ylabel('Skewness', fontsize=20)
    r2text = r"${}^{{2}} = {:.{p}f} $".format('r',r_value**2 , p=2)
    plt.text(0.1, 0.62 , r'$p < 0.00001$', fontsize=14, horizontalalignment='center', \
                verticalalignment='center',transform = ax.transAxes)
    plt.text(0.1, 0.67, r2text,  fontsize=14, horizontalalignment='center', \
                verticalalignment='center',transform = ax.transAxes)
    fig.savefig(git_path + 'figs/RSG.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

plotRSG()
