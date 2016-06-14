import sys, os, copy
from Bio import SeqIO

mydir = os.path.expanduser("~/github/Task2/LTDE")


def annotateGATK(species, isopen=False, CUTOFF=25.0, start=3, stop=4, name=8, strand=6):
    GBK_path = mydir + '/data/2015_SoilGenomes_Annotate/' + str(species) + '/G-Chr1.gbk'
    print GBK_path
    GBK = SeqIO.parse(open(GBK_path,"r"), "genbank")
    for gb_record in GBK:
    	contig = gb_record.name
    	contig_seq = gb_record.seq
        gb_feature = gb_record.features
        for gb in gb_feature:
            if gb.type == "CDS" and "product" in gb.qualifiers:
                print gb.qualifiers
                print gb.qualifiers['product']
    #print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
    #print repr(gb_record.seq)

strains = ['KBS0703', 'KBS0705', 'KBS0706', 'KBS0710', 'KBS0711', 'KBS0713', \
    'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802']

for x in strains:
    annotateGATK(x)
