from __future__ import division
import FlowCytometryTools
from FlowCytometryTools import test_data_dir, test_data_file, FCMeasurement, FCPlate, ThresholdGate, PolyGate
import os, re, datetime, signal, collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
    if analysis == 'RSG':
        for strain in strains:
            for transfer in transfers:
                path = git_path + 'figs/' + analysis + '/' + strain + '/' + transfer + '/'
                if not os.path.exists(path):
                    os.makedirs(path)
    else:
        for strain in strains:
            for transfer in transfers:
                for rep in reps:
                    rep_path = git_path + 'figs/' + analysis + '/' + strain + '/' + transfer + '/' + rep + '/'
                    if not os.path.exists(rep_path):
                        os.makedirs(rep_path)

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

def get_KDE_RSG():
    makeFigFolders(analysis = 'RSG')
    strains = ['B', 'C', 'D', 'F', 'P', 'J']
    transfers = ['D1', 'D10', 'D100']
    reps = ['1', '2', '3', '4', '5']
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

                    b = datetime.datetime.now()
                    c = b - a
                    print str(c.seconds) + " seconds"

                except TimeoutException:
                    continue # continue the for loop if function takes more than x seconds
                else:
                    # Reset the alarm
                    signal.alarm(0)

def getDistRSG():
    makeFigFolders(analysis = 'RSG')
    strains = ['B', 'C', 'D', 'F', 'P', 'J']
    transfers = ['D1', 'D10', 'D100']
    reps = ['1', '2', '3', '4', '5']
    dilution_dict = {'D1': str(0), 'D10': str(1), 'D100':str(2)}
    OUT1 = open(git_path +'data/FlowSummary.txt', 'w+')
    print>> OUT1, 'Strain', 'Transfer', 'Rep', 'Mean', 'Std', 'Median', 'Skew', 'cells_per_ml'
    dilution_df = pd.read_excel(git_path +'data/Flow_dilutions/Flow_dilutions_D100_test.xlsx')
    s = dilution_df['Sample'].str.split(' ').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1) # to line up with df's index
    s.name = 'Sample' # needs a name to join
    dilution_df_split = dilution_df.join(s.apply(lambda x: pd.Series(x.split('-'))))

    dilution_df_split.rename(columns={0: 'Transfer', 1: 'Shaker', 2: 'Position', \
    3: 'Dilution', 4: 'Strain', 5:'Replicate'},inplace=True)
    for strain in strains:
        for transfer in transfers:
            for rep in reps:
                path = data_path + strain + '/' + transfer + '/' + rep + '/'
                files = os.listdir(path)
                remove = ['misc', '.DS_Store']
                files1 = [i for i in files if i not in remove]
                files_split = [re.split(r'[._]',x)[5] for x in files1]
                if 'I' in files_split:
                    dilution_df_split_line = dilution_df_split[ \
                        (dilution_df_split.Transfer == dilution_dict[transfer]) & \
                        (dilution_df_split.Strain == strain) & \
                         (dilution_df_split.Replicate == rep)]
                    Flow_dilution = dilution_df_split_line.Flow_dilution

                    path = data_path + strain + '/' + transfer + '/' + rep + '/'
                    plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
                    plate = plate.dropna()
                    plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
                        'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
                    DAPI_threshold = getDAPIgate(plate)
                    DAPI_gate = ThresholdGate(DAPI_threshold, 'Pacific Blue-A', region='above')
                    cell_gate = ThresholdGate(250000, 'SSC-H', region='below')
                    bead_gate = ThresholdGate(250000, 'SSC-H', region='above')
                    intersection_gate = DAPI_gate & cell_gate
                    # cell count
                    DAPI_gated_A9 = plate['A9'].gate(intersection_gate)
                    cells = len(DAPI_gated_A9.data[['SSC-H']].values)
                    # bead count
                    beads_gated_A9 = plate['A9'].gate(bead_gate)
                    beads = len(beads_gated_A9.data[['SSC-H']].values)
                    # RSG
                    RSG = DAPI_gated_A9.data[['FITC-A']].values
                    # (total cells / total beads) = (counted cells / counted beads)
                    # total cells = (counted cells / counted beads) * total beads
                    # total beads = 10 uL of 1x10^8 beads/mL = 1*10^8 * 0.01 = 1*10^6 beads
                    total_cells = (cells / beads) * 1000000
                    # Flow_dilution = fold dilution (ex. 10 = 10-fold dilution)
                    cells_per_ml = total_cells * Flow_dilution
                    cells_per_ml = cells_per_ml.values[0]
                    print strain, transfer, rep, str(cells_per_ml)

                    print>> OUT1, strain, transfer[1:], rep, str(np.mean(RSG)), str(np.std(RSG)), \
                        str(np.median(RSG)), str(stats.skew(RSG)[0]), str(cells_per_ml)

    OUT1.close()


def testCounts():
    strains = ['B', 'C', 'D', 'F', 'P', 'J']
    transfers = ['D1', 'D10', 'D100']
    reps = ['1', '2', '3', '4', '5']
    samples = ['A8', 'A9']
    makeFigFolders('beads')
    for strain in strains:
        for transfer in transfers:
            for rep in reps:
                print strain, transfer, rep
                path = data_path + strain + '/' + transfer + '/' + rep + '/'
                files = os.listdir(path)
                remove = ['misc', '.DS_Store']
                files1 = [i for i in files if i not in remove]
                files_split = [re.split(r'[._]',x)[5] for x in files1]
                if 'I' in files_split:
                    for sample in samples:
                        plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
                        plate = plate.dropna()
                        plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
                            'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
                        threshold = getDAPIgate(plate)
                        gate = ThresholdGate(threshold, 'Pacific Blue-A', region='above')
                        gated_sample = plate[sample].gate(gate)
                        FSC_A = gated_sample.data[['FSC-A']].values
                        SSC_A = gated_sample.data[['SSC-A']].values
                        PI = gated_sample.data[['PI (YG)-A']].values
                        FITC = gated_sample.data[['FITC-A']].values
                        FSC_H = gated_sample.data[['FSC-H']].values
                        SSC_H = gated_sample.data[['SSC-H']].values
                        FSC_PMT_A = gated_sample.data[['FSC PMT-A']].values
                        FSC_PMT_H = gated_sample.data[['FSC PMT-H']].values
                        Pacific_blue = gated_sample.data[['Pacific Blue-A']].values
                        APC = gated_sample.data[['APC-A']].values

                        fig = plt.figure()
                        ax = fig.add_subplot(111, projection='3d')
                        ax.scatter(FITC, Pacific_blue, APC , zdir='z', alpha = 0.2)
                        ax.set_xlabel('FITC')
                        ax.set_ylabel('Pacific blue')
                        ax.set_zlabel('APC')


                        name = sample + '_FITC_APC_PacificBlue_' + strain + '_' + transfer + '_' + rep
                        plt.savefig(git_path + 'figs/beads'  + '/' + strain + '/' \
                            + transfer + '/' + rep + '/' + name + '.png')
                        plt.close()

                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.scatter(APC, FITC, alpha = 0.2)
                        ax.set_xlabel('APC')
                        ax.set_ylabel('FITC')

                        name = sample +'_APC_FITC_' + strain + '_' + transfer + '_' + rep
                        plt.savefig(git_path + 'figs/beads'  + '/' + strain + '/' \
                            + transfer + '/' + rep + '/' + name + '.png')
                        plt.close()

                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.scatter(Pacific_blue, FITC, alpha = 0.2)
                        ax.set_xlabel('Pacific blue')
                        ax.set_ylabel('FITC')

                        name = sample + '_PacificBlue_FITC_' + strain + '_' + transfer + '_' + rep
                        plt.savefig(git_path + 'figs/beads'  + '/' + strain + '/' \
                            + transfer + '/' + rep + '/' + name + '.png')
                        plt.close()

                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.scatter(FSC_A, SSC_A, alpha = 0.2)
                        ax.set_xlabel('FSC-A')
                        ax.set_ylabel('SSC-A')

                        name = sample + '_FSC_A_SSC_A_' + strain + '_' + transfer + '_' + rep
                        plt.savefig(git_path + 'figs/beads'  + '/' + strain + '/' \
                            + transfer + '/' + rep + '/' + name + '.png')
                        plt.close()

                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.scatter(FSC_H, SSC_H, alpha = 0.2)
                        ax.set_xlabel('FSC-H')
                        ax.set_ylabel('SSC-H')
                        #plt.axhline(y = 250000, linewidth=2, color='darkgrey', ls='--')

                        name = sample + '_FSC_H_SSC_H_' + strain + '_' + transfer + '_' + rep
                        plt.savefig(git_path + 'figs/beads'  + '/' + strain + '/' \
                            + transfer + '/' + rep + '/' + name + '.png')
                        plt.close()

                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.scatter(FSC_PMT_A, FSC_PMT_H, alpha = 0.2)
                        ax.set_xlabel('FSC PMT-A')
                        ax.set_ylabel('FSC PMT-H')

                        name = sample + '_FSC_PMT_A_FSC_PMT_H_' + strain + '_' + transfer + '_' + rep
                        plt.savefig(git_path + 'figs/beads'  + '/' + strain + '/' \
                            + transfer + '/' + rep + '/' + name + '.png')
                        plt.close()



def plotRSG():
    strain_dict = {'A':['cyan', 'Arthrobacter'], 'B':['lightblue', 'Bacillus'], \
            'C':['lightgreen', 'Caulobacter'], 'D':['darkred', 'Deinococcus'], 'F':['orange', 'Pedobacter'], \
            'J':['indigo', 'Janthinobacterium'], 'P':['darkgreen', 'Pseudomonas']}
    IN = pd.read_csv(git_path +'data/bandwidth.txt', sep = ' ')
    x = np.log10(IN.Transfer.values)
    y = IN.Mean.values
    y_std = [ (y_i - np.mean(y)) / np.std(y) for y_i in y ]
    IN['MeanStd'] = pd.Series(y_std, index=IN.index)
    IN = IN.sort_values(by=["Transfer","MeanStd"], ascending=[True, False])
    names_count = collections.Counter(IN.Strain.tolist())
    names_count_sorted  = sorted(names_count.items(), key=itemgetter(0))

    groups = IN.groupby('Strain')
    # Plot
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in groups:
        ax.plot(np.log10(group.Transfer), group.MeanStd, marker='o', alpha = 0.4 , \
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
    plt.tick_params(axis='both', which='major', labelsize=10)
    #plt.ylim([min(y)-0.1,max(y)+0.1])
    #plt.xlim([min(x)-0.001,max(x)+0.001])
    plt.xlabel('Days starved, log10', fontsize=20)
    plt.ylabel('Standardized mean', fontsize=20)
    r2text = r"${}^{{2}} = {:.{p}f} $".format('r',r_value**2 , p=2)
    #plt.text(0.1, 0.62 , r'$p < 0.00001$', fontsize=14, horizontalalignment='center', \
    #            verticalalignment='center',transform = ax.transAxes)
    #plt.text(0.1, 0.67, r2text,  fontsize=14, horizontalalignment='center', \
    #             verticalalignment='center',transform = ax.transAxes)
    fig.savefig(git_path + 'figs/RSG.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def plotStrainsRSG(moment = 'Mean'):
    # accepts 'Mean', 'Std', and 'Skew' as arguments for moment
    strain_dict = {'A':['cyan', 'Arthrobacter'], 'B':['lightblue', 'Bacillus'], \
            'C':['lightgreen', 'Caulobacter'], 'D':['darkred', 'Deinococcus'], 'F':['orange', 'Pedobacter'], \
            'J':['indigo', 'Janthinobacterium'], 'P':['darkgreen', 'Pseudomonas']}
    IN = pd.read_csv(git_path +'data/FlowSummary.txt', sep = ' ')
    fig = plt.figure()
    count = 0
    for strain in strain_dict.keys():
        if strain  == 'A':
            continue
        else:
            IN_strain = IN.loc[IN['Strain'] == strain]
            x = np.log10(IN_strain.Transfer.values)
            if moment == 'Mean':
                y = IN_strain.Mean.values
                #y1 = FlowCytometryTools.core.transforms.hlog_inv(y1)
            elif moment == 'Std':
                y = IN_strain.Std.values
            elif moment == 'Skew':
                y = IN_strain.Skew.values

            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            #print strain_dict[strain][1]
            #print "slope = " + str(slope)
            #print "r2 = " + str(r_value**2)
            #print "p = " + str(p_value)
            predict_y = intercept + slope * x
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
            ax = fig.add_subplot(2,3,count + 1)
            plt.plot(x, predict_y, 'k-')
            ax.scatter(x, y)
            ax.set_title(strain_dict[strain][1])
            ax.tick_params(labelsize=6)
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
            count += 1
    fig.text(0.50, 0.001, 'Starvation time, ' + r'$log_{10}$', ha='center', va='center', fontsize=16)
    fig.text(0.02, 0.5,  moment + ' reductase activity', ha='center', va='center', rotation='vertical', fontsize=16)
    fig.savefig(git_path + 'figs/RSG_Strains' + moment + '.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plotStrainsCells():
    strain_dict = {'A':['cyan', 'Arthrobacter'], 'B':['lightblue', 'Bacillus'], \
            'C':['lightgreen', 'Caulobacter'], 'D':['darkred', 'Deinococcus'], 'F':['orange', 'Pedobacter'], \
            'J':['indigo', 'Janthinobacterium'], 'P':['darkgreen', 'Pseudomonas']}
    IN = pd.read_csv(git_path +'data/FlowSummary.txt', sep = ' ')
    fig = plt.figure()
    count = 0
    for strain in strain_dict.keys():
        if strain  == 'A':
            continue
        else:
            IN_strain = IN.loc[IN['Strain'] == strain]
            x = np.log10(IN_strain.Transfer.values)
            y = np.log10(IN_strain.cells_per_ml.values)
            print strain
            print x
            print y
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            #print strain_dict[strain][1]
            #print "slope = " + str(slope)
            #print "r2 = " + str(r_value**2)
            #print "p = " + str(p_value)
            predict_y = intercept + slope * x
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
            ax = fig.add_subplot(2,3,count + 1)
            plt.plot(x, predict_y, 'k-')
            ax.scatter(x, y)
            ax.set_title(strain_dict[strain][1])
            ax.tick_params(labelsize=6)
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
            count += 1
    fig.text(0.50, 0.001, 'Starvation time, ' + r'$log_{10}$', ha='center', va='center', fontsize=16)
    fig.text(0.02, 0.5, 'Cells per mL, ' + r'$log_{10}$', ha='center', va='center', rotation='vertical', fontsize=16)
    fig.savefig(git_path + 'figs/Cells_Strains.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def plotStrainsCellRSG(moment = 'Mean'):
    # accepts 'Mean', 'Std', and 'Skew' as arguments for moment
    strain_dict = {'A':['cyan', 'Arthrobacter'], 'B':['lightblue', 'Bacillus'], \
            'C':['lightgreen', 'Caulobacter'], 'D':['darkred', 'Deinococcus'], 'F':['orange', 'Pedobacter'], \
            'J':['indigo', 'Janthinobacterium'], 'P':['darkgreen', 'Pseudomonas']}
    IN = pd.read_csv(git_path +'data/FlowSummary.txt', sep = ' ')
    fig = plt.figure()
    count = 0
    for strain in strain_dict.keys():
        if strain  == 'A':
            continue
        else:
            IN_strain = IN.loc[IN['Strain'] == strain]
            x = np.log10(IN_strain.Transfer.values)
            if moment == 'Mean':
                y1 = IN_strain.Mean.values
                #y1 = FlowCytometryTools.core.transforms.hlog_inv(y1)
            elif moment == 'Std':
                y1 = IN_strain.Std.values
            elif moment == 'Skew':
                y1 = IN_strain.Skew.values
            y2 = np.log10(IN_strain.cells_per_ml.values)

            y = y1 / y2
            print strain

            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            #print strain_dict[strain][1]
            #print "slope = " + str(slope)
            #print "r2 = " + str(r_value**2)
            #print "p = " + str(p_value)
            predict_y = intercept + slope * x
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
            ax = fig.add_subplot(2,3,count + 1)
            plt.plot(x, predict_y, 'k-')
            ax.scatter(x, y)
            ax.set_title(strain_dict[strain][1])
            ax.tick_params(labelsize=6)
            plt.setp(ax.get_xticklabels()[::2], visible=False)
            plt.setp(ax.get_yticklabels()[::2], visible=False)
            count += 1
    fig.text(0.50, 0.001, 'Starvation time, ' + r'$log_{10}$', ha='center', va='center', fontsize=16)
    fig.text(0.02, 0.5,  moment + ' reductase activity/ cell/ mL', ha='center', va='center', rotation='vertical', fontsize=16)
    fig.savefig(git_path + 'figs/RSG_Cell_Strains' + moment + '.png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

#getDistRSG()
plotStrainsRSG()
plotStrainsRSG(moment = 'Skew')

plotStrainsCells()
plotStrainsCellRSG()
