from __future__ import division
import sys, os, copy, re
import pandas as pd

class annotate:

    '''
    Feed the class the species name
    '''

    def __init__ (self, GFF, FFN):
        self.GFF = GFF
        self.FFN = FFN
        self.start =3
        self.stop = 4
        self.name = 8
        self.strand = 6
        self.CUTOFF=25.0
        self.isopen=False
        self.codon_dict = {
    "TTT":"F", "TCT":"S", "TAT":"Y", "TGT":"C",
    "TTC":"F", "TCC":"S", "TAC":"Y", "TGC":"C",
    "TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*",
    "TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W",

    "CTT":"L", "CCT":"P", "CAT":"H", "CGT":"R",
    "CTC":"L", "CCC":"P", "CAC":"H", "CGC":"R",
    "CTA":"L", "CCA":"P", "CAA":"Q", "CGA":"R",
    "CTG":"L", "CCG":"P", "CAG":"Q", "CGG":"R",

    "ATT":"I", "ACT":"T", "AAT":"N", "AGT":"S",
    "ATC":"I", "ACC":"T", "AAC":"N", "AGC":"S",
    "ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R",
    "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R",

    "GTT":"V", "GCT":"A", "GAT":"D", "GGT":"G",
    "GTC":"V", "GCC":"A", "GAC":"D", "GGC":"G",
    "GTA":"V", "GCA":"A", "GAA":"E", "GGA":"G",
    "GTG":"V", "GCG":"A", "GAG":"E", "GGG":"G"
    }

    class gene:
        def __init__ (self, name, start, stop, strand):
            self.name=name
            self.start=start
            self.stop=stop
            self.strand=strand

    def is_gene(self, line):
        try:
            if line[2]=="CDS":
                return True
            else:
                return False
        except:
            return False

    def is_fourfold(self, test):
        fourfold_first_two = ['CG', 'GC', 'GG', 'CT', 'CC', 'TC', 'AC', 'GT']
        if test[:-1] in fourfold_first_two:
            return 'Y'
        else:
            return 'N'


    def rcomp(self, seq):
        s=[]
        rc={'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        rseq=[]
        for s in seq:
            rseq.insert(0, rc[s])
        return ''.join(rseq)

    def getmut(self, start, bp, seq, mut):
        A=int(bp)-(int(bp)-int(start))%3
        B=int(bp)
        return (self.codon_dict[seq[A:A+3]]+str((bp-start)/3+1)+\
            self.codon_dict[seq[A:B]+mut+seq[B+1:A+3]], seq[A:A+3], seq[A:B]+mut+seq[B+1:A+3] )

    def getrev(self, start, bp, seq, mut):
        bp=int(bp)+1
        start=int(start)+1

        A=int(bp)+abs(int(start)-int(bp))%3
        B=int(bp)
        return (self.codon_dict[self.rcomp(seq[A-3:A])]+str((start-bp)/3+1)+self.codon_dict[self.rcomp(seq[A-3:B-1]+mut+seq[B:A])], \
            self.rcomp(seq[A-3:A]), self.rcomp(seq[A-3:B-1]+mut+seq[B:A]))

    def getGenesAndSizes(self, GFF):
        gene_sizes={}
        genes=[]
        for line in GFF:
            line=line.split('\t')
            if self.is_gene(line):
                for x in line[self.name].split(';'):
                    if x.split('=')[0]=="gene":
                        genes.append(annotate.gene(x.split('=')[1], int(line[self.start]), int(line[self.stop]), line[self.strand]) )
                        gene_sizes[x.split('=')[1]]=(int(line[self.start])-int(line[self.stop]))
        return (gene_sizes, genes)

    def getSeqs(self, FFN):
        seq=[]
        for line in FFN:
            if line[0]!='>':
                line=line.strip('\n')
                seq.append(line)
        seq=' '+''.join(seq)
        return seq

    def annotateGATK(self, IN, OUT):
        IN = open(IN)
        OUT = open(OUT, 'wr')
        GFF = open(self.GFF)
        FFN = open(self.FFN)
        get_tuple = self.getGenesAndSizes(GFF)
        gene_sizes = get_tuple[0]
        genes = get_tuple[1]
        this_gene=genes.pop(0)
        seq = self.getSeqs(FFN)

        this_gene_copy = copy.copy(this_gene)
        genes_copy = copy.copy(genes)
        seq_copy = copy.copy(seq)
        for line in IN:
            fields=line.strip().split('\t')
            if fields[0]=='CHROM':
                fields.insert( 5, 'GENE')
                fields.insert( 6, 'CODING')
                fields.insert( 7, 'FOURFOLD')
                print>>OUT, fields[0], fields[1], fields[3], fields[4], \
                    fields[5], fields[6], fields[7], fields[8], fields[9]
                continue
            ref=fields[3]
            alt=fields[4]
            pos=int(fields[1])
            if len(alt)==1:
                while ((pos>max(this_gene_copy.start, this_gene_copy.stop ))  and len(genes_copy)  >0 ):
                    this_gene_copy=genes_copy.pop(0)
                if this_gene_copy.strand=='+':
                    if (pos>this_gene_copy.start):
                        getmut=self.getmut(this_gene_copy.start, pos, seq, alt)
                        kind = getmut[0]
                        hit=(this_gene_copy.name+":"+kind)
                        if kind[0]==kind[-1]:
                            kind='S'
                            fourfold = self.is_fourfold(getmut[1])
                        else:
                            kind='N'
                            fourfold = 'N'
                    else:
                        hit='NC'
                        kind='NC'
                        fourfold = 'NC'
                elif this_gene_copy.strand=='-':
                    if (pos>this_gene_copy.start):
                        getrev = self.getrev(this_gene_copy.stop, pos, seq, alt)
                        kind = getrev[0]
                        hit=(this_gene_copy.name+":"+kind)
                        if kind[0]==kind[-1]:
                            kind='S'
                            fourfold = self.is_fourfold(getrev[1])
                        else:
                            kind='N'
                            fourfold = 'N'
                    else:
                        hit='NC'
                        kind='NC'
                        fourfold = 'NC'
                fields.insert(5, hit)
                fields.insert(6, kind)
                fields.insert(7, fourfold)
                print>>OUT, fields[0], fields[1], fields[3], fields[4], \
                    fields[5], fields[6], fields[7], fields[8],  fields[9]

        OUT.close()

    def annotateMAPGD(self, IN, OUT):
        IN = open(IN)
        OUT = open(OUT, 'w')
        GFF = open(self.GFF)
        FFN = open(self.FFN)
        get_tuple = self.getGenesAndSizes(GFF)
        gene_sizes = get_tuple[0]
        genes = get_tuple[1]
        this_gene=genes.pop(0)
        seq = self.getSeqs(FFN)

        this_gene_copy = copy.copy(this_gene)
        genes_copy = copy.copy(genes)
        seq_copy = copy.copy(seq)

        for line in IN:
            if self.isopen:
                if line[0]=='@':
                    print line,
                    self.isopen=False
                fields=line.split('\t')
                ChiPoly=[]
                ChiFixed=[]
                Freq=[]
                for X in fields[6:]:
                    num=X.split('/')
                    try:
                        ChiPoly.append(float(num[2]))
                        ChiFixed.append(float(num[3]))
                        P=float(num[0])
                        cov=int(num[1])
                        if ChiPoly[-1]<self.CUTOFF:
                            P=round(P)
                        if cov>0:
                            Freq.append(P)
                        else:
                            Freq.append('N')
                    except:
                        ChiPoly.append(0)
                        ChiFixed.append(0)
                        Freq.append('N')
                if fields[0][0] == '@':
                    continue
                ref=fields[2]
                alt=fields[3]
                pos=int(fields[1])
                if len(alt)==1:
                    while (pos>max(this_gene.start, this_gene.stop ) ):
                        this_gene=genes.pop(0)

                    if this_gene.strand=='+':
                        if (pos>this_gene.start):
                            getmut=self.getmut(this_gene.start, pos, seq, alt)
                            kind = getmut[0]
                            hit=(this_gene.name+":"+kind)
                            if kind[0]==kind[-1]:
                                kind='S'
                                fourfold = self.is_fourfold(getmut[1])
                            else:
                                kind='N'
                                fourfold = 'N'
                        else:
                            hit='NC'
                            kind='NC'
                            fourfold = 'NC'
                    elif this_gene.strand=='-':
                        if (pos>this_gene.start):
                            getrev=self.getrev(this_gene.stop, pos, seq, alt)
                            kind = getrev[0]
                            hit=(this_gene.name+":"+kind)
                            if kind[0]==kind[-1]:
                                kind='S'
                                fourfold = self.is_fourfold(getrev[1])
                            else:
                                kind='N'
                                fourfold = 'N'
                        else:
                            hit='NC'
                            kind='NC'
                            fourfold = 'NC'
                if max(ChiPoly)>self.CUTOFF or max(ChiFixed)>self.CUTOFF:
                    print '\t'.join(map(str, fields[:6]))+'\t'+hit+'\t'+kind+'\t'+ fourfold +'\t'+'\t'.join(map(str, Freq))
                    print>> OUT, '\t'.join(map(str, fields[:6]))+'\t'+hit+'\t'+kind +'\t'+ fourfold +'\t'+'\t'.join(map(str, Freq))
            else:
                print line,
                fields=line.split('\t')
                if fields[0]=="@SCFNAME       ":
                    self.isopen=True
        OUT.close

class formatSNPs:

    def __init__(self, strain):
        #self.pathIN = pathIN
        #self.pathOUT = pathOUT
        self.strain = strain
        self.mydir = os.path.expanduser("~/github/Task2/LTDE")

    def getSNPTable(self):
        dfs = []
        samples = []
        pathIN = self.mydir + '/data/GATK/annotate/' + self.strain
        pathOUT = self.mydir + '/data/GATK/merged/all/'  + self.strain
        for i in os.listdir(pathIN):
            if i.endswith(".txt"):
                sample_name = i.split('.')[0]
                sample = re.split(r'[-_]+', sample_name)[2]
                samples.append(sample)
                names = ['Scaffold', 'Pos', 'Ref', 'Alt', 'Gene', \
                    'Coding', 'Fourfold', 'Qual', 'AC']
                IN = pd.read_csv(pathIN + '/' + i, delimiter = ' ', header = None)
                IN.columns = names
                # remove original header and AC column
                IN = IN.iloc[1:,:-1]
                dfs.append(IN)
            else:
                continue
        merged_df = reduce(lambda left,right: pd.merge(left,right,on=['Scaffold', 'Pos','Ref','Gene'], how='outer'), dfs)
        toRename = ['Alt', 'Coding', 'Fourfold','Qual']
        new_columns = ['Scaffold', 'Pos', 'Ref', 'Gene']
        to_move = ''
        for x, y in enumerate(samples):
            for z in toRename:
                out = z + '_' + y
                if x ==0 and z == 'Alt':
                    new_columns.insert(3,out)
                    to_move = out
                else:
                    new_columns.append(out)
        merged_df.columns = new_columns
        to_move_values = merged_df.loc[:,to_move]
        merged_df.drop(to_move, axis=1, inplace=True)
        merged_df.insert(4, to_move, to_move_values)
        merged_df.to_csv(pathOUT +'_GATK.txt' ,sep='\t', \
            index = False)

    def filterSNPTable(self, coding_to_csv=False):
        pathIN = self.mydir + '/data/GATK/merged/all/'  + self.strain
        IN = pd.read_csv(pathIN + '_GATK.txt', delimiter = '\t', header = 'infer')
        IN_coding = IN[IN['Gene'] != 'NC']
        if coding_to_csv == True:
            IN_coding.to_csv(self.mydir + '/data/GATK/merged/coding/' + self.strain +'_GATK_C.txt' ,sep='\t', \
                index = False)
        NAN_counts = IN_coding.isnull().sum(axis=1)
        toKeep = NAN_counts[NAN_counts != 0].index
        toKeep_df = IN_coding.ix[toKeep]
        toKeep_df.to_csv(self.mydir + '/data/GATK/merged/potential_coding_mutations/' + self.strain +'_GATK_C_PCM.txt' ,sep='\t', \
            index = False)


def annotateStrains():
    mydir = os.path.expanduser("~/github/Task2/LTDE")
    strains = ['KBS0703', 'KBS0705', 'KBS0706', 'KBS0710', 'KBS0711', 'KBS0713', \
        'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802', 'KBS0812']
    #strains = ['KBS0711']
    for strain in strains:
        content_list = []
        for content in os.listdir(mydir + '/data/GATK/raw/' + strain): # "." means current directory
            content_list.append(content)
        if strain == 'KBS0812':
            GFF = mydir + '/data/Bacillus_test/AL009126.3.gff'
            FFN = mydir + '/data/Bacillus_test/AL009126.fa'
        else:
            GFF = mydir + '/data/2015_SoilGenomes_Annotate/' + strain + '/G-Chr1.gff'
            FFN = mydir + '/data/2015_SoilGenomes_Annotate/' + strain + '/G-Chr1.fna'
        mapgd=mydir + '/data/mapgd/raw/' + str(strain) + '_merged.pol'
        mapgdOUT =  mydir+'/data/mapgd/annotate/' + str(strain) +  '_merged_annotate.txt'
        annotate(GFF, FFN).annotateMAPGD(mapgd, mapgdOUT)
        for x in content_list:
            IN = mydir + '/data/GATK/raw/' +  strain + '/' + x
            OUT =  mydir + '/data/GATK/annotate/' + str(strain) + '/' + x.split('.')[0] + '_annotate.txt'
            annotate(GFF, FFN).annotateGATK(IN, OUT)

        formatSNPs(strain).getSNPTable()
        formatSNPs(strain).filterSNPTable(coding_to_csv=True)

annotateStrains()
