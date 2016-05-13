from __future__ import division
import os, argparse, random, math, re, collections
from collections import Counter
from itertools import product
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
import pandas as pd
import numpy as np

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'


def importAAtable(AAfile):
    table = {}
    aa_list = [[item.strip() for item in line.rstrip('\r\n').split('= ', 1)[-1]] \
        for line in open(AAfile)]
    aa_list_zip = map(list, zip(*aa_list))
    for x in aa_list_zip:
         x[2:5] = [''.join(x[2:5])]
         table[str(x[2])] = str(x[0])
    return table

test = mydir + 'codon.txt'
AAtable = importAAtable(test)

# First, parse.pol file

def parsePOL(pol):
    for line in open(pol):
        if line.startswith('@'):
            test = line.split()
            print test

pol = mydir + 'mapgd/GSF1018_merged_sorted_NOdup_sorted.pol'

polIN = pd.read_csv(filepath_or_buffer= pol, sep ='\t', skiprows = 2, header = None)
# drop NaN column
polIN = polIN.drop(polIN.columns[[7]], axis=1)
# split last column
#polIN.join(s.apply(lambda x: Series(x.split(':'))))
polSplit = polIN[polIN.columns[-1]].str[1:-1].str.split('/', return_type='frame').astype(float)
polIN = polIN.drop(polIN.columns[[6]], axis=1)
frames = [polIN, polSplit]
result = pd.concat(frames, axis=1)
# add headers
###NOTEs: first item in name contains locus tag
result.columns = ['Name', 'Position', 'Major', 'Minor', 'Coverage', 'Error', 'MLE', 'LLRT', 'x1', 'x2']
#### import gbk

contigs = [ int(re.search(r'\d+', x).group()) for x in result['Name'].tolist()]

contigsPositions = [contigs, result['Position'].tolist(), \
        result['Major'].tolist(), result['Minor'].tolist()]
#print len(set(contigsPositions[0]))
#print contigsPositions[0]

def returnCodon(start, end, x, feature_seq, contigsPositions):
    '''
    codonPosition = 0, 1, or 2
    feature_seq = Sequence of CDS from biopython
    contigsPositions = nested list containing = [[codons], [positions],
        [major allele], [minor allele]]
    '''
    #print feature_seq[contigsPositions[1][x]]
    length_seq = end - start
    position = contigsPositions[1][x] - 1 - start
    print contigsPositions[1][x], length_seq, position
    diff =  length_seq - contigsPositions[1][x] -1
    codonPosition =  position % 3
    if codonPosition == 0:
        codon = feature_seq[position:position+3]
        minor_codon = contigsPositions[2][x] + \
            (feature_seq[position+1:position+3])
        major_codon = contigsPositions[3][x] + \
            (feature_seq[position+1:position+3])
    if codonPosition == 1:
        codon = feature_seq[position-1:position+2]
        minor_codon = Seq(str(feature_seq[position-1] + \
            contigsPositions[2][x] + \
            (feature_seq[position+2])), generic_dna)
        major_codon = Seq(str(feature_seq[position-1] + \
            contigsPositions[3][x] + \
            (feature_seq[position+2])), generic_dna)
    if codonPosition == 2:
        codon = feature_seq[position-2:position+1]
        minor_codon = (feature_seq[position-2:x]) + \
            contigsPositions[2][x]
        major_codon = (feature_seq[position-2:x]) + \
            contigsPositions[3][x]
    returnTuple = (codon, major_codon, minor_codon, codonPosition)
    return returnTuple

def majorMinorCodon(codon, major_codon, minor_codon):
    if major_codon == codon:
        Reference = 'MJ'
        Alternative = 'MN'
    elif minor_codon == codon:
        Reference = 'MN'
        Alternative = 'MJ'
    else:
        Reference = 'NA'
        Alternative = 'NA'
    #if len(codon) != 3 and \
    #    len(major_codon) != 3 and \
    #    len(minor_codon) != 3:
    #    pass
    returnTuple = (Reference, Alternative)
    return returnTuple

def GBK_parse(gbk, contigsPositions):
    count = 0
    BioPy_gbk = SeqIO.parse(gbk, "genbank")
    gb_dict = {}
    for record in BioPy_gbk:
        contigID = int(re.search(r'\d+', str(record.id.strip())).group())
        if contigID not in contigsPositions[0]:
            continue
        else:
            listOccur = [i for i, x in enumerate(contigsPositions[0]) if x == contigID]
            for x in listOccur:
                position = contigsPositions[1][x]
                types = 0
                for feature in record.features:
                    mystart = feature.location._start.position
                    site = feature.location
                    myend = feature.location._end.position
                    #print feature.type
                    if position in range(mystart, myend) and count not in gb_dict \
                        and (feature.type == 'CDS' ) :
                        feature_seq = feature.extract(record.seq)
                        product = feature.qualifiers['product']
                        returnTuple = returnCodon(mystart, myend, x, \
                            feature_seq, contigsPositions)

                        codon = returnTuple[0]
                        major_codon = returnTuple[1]
                        minor_codon = returnTuple[2]
                        codonPosition =returnTuple[3]
                        #majorMinor = majorMinorCodon(codon, major_codon, minor_codon)
                        Reference = feature_seq[x]
                        snp_0_tl = str(codon.translate())
                        snp_MJ_tl = str(major_codon.translate())
                        snp_MN_tl = str(minor_codon.translate())
                        gb_dict[count] = snp_0_tl, \
                            snp_MJ_tl, snp_MN_tl, codonPosition, product
                        #gb_dict[count] = Reference, Alternative, \
                        #    snp_MJ_tl, snp_MN_tl, codonPosition
                        count +=1
    print gb_dict
    returnData = pd.DataFrame.from_dict(gb_dict, orient='index')
    returnData.columns = ['Reference', 'MajorAlleleTrans', \
            'MinorAlleleTrans', 'codonPosition', 'product']
    #returnData.columns = ['Reference', 'Alternative', 'MajorAlleleTrans', \
    #        'MinorAlleleTrans', 'codonPosition']
    return returnData

gbk = mydir + 'reference_genomes/X_03292016.gbk'


MAPGDannotate = GBK_parse(gbk, contigsPositions)
print MAPGDannotate
filesToMerge = [result, MAPGDannotate]
mergePandas = pd.concat(filesToMerge, axis=1)

splittt = pol.split('/')[-1]
splitttt = splittt.split('.')[-1]

outAnnotate = mydir + 'mapgd/' + splitttt + '_annotate.txt'
#path = mydir +
mergePandas.to_csv(outAnnotate)
#print mergePandas
