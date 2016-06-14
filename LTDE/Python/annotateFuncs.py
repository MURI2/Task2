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
        return self.codon_dict[seq[A:A+3]]+str((bp-start)/3+1)+self.codon_dict[seq[A:B]+mut+seq[B+1:A+3]]

    def getrev(self, start, bp, seq, mut):
        bp=int(bp)+1
        start=int(start)+1

        A=int(bp)+abs(int(start)-int(bp))%3
        B=int(bp)
        return self.codon_dict[self.rcomp(seq[A-3:A])]+str((start-bp)/3+1)+self.codon_dict[self.rcomp(seq[A-3:B-1]+mut+seq[B:A])]

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
                print>>OUT, fields[0], fields[1], fields[3], fields[4], \
                    fields[5], fields[6], fields[7], fields[8]
                continue
            ref=fields[3]
            alt=fields[4]
            pos=int(fields[1])
            if len(alt)==1:
                while ((pos>max(this_gene_copy.start, this_gene_copy.stop ))  and len(genes_copy)  >0 ):
                    this_gene_copy=genes_copy.pop(0)
                if this_gene_copy.strand=='+':
                    if (pos>this_gene_copy.start):
                        kind=self.getmut(this_gene_copy.start, pos, seq, alt)
                        hit=(this_gene_copy.name+":"+kind)
                        if kind[0]==kind[-1]:
                            kind='S'
                        else:
                            kind='N'
                    else:
                        hit='NC'
                        kind='NC'
                elif this_gene_copy.strand=='-':
                    if (pos>this_gene_copy.start):
                        kind = self.getrev(this_gene_copy.stop, pos, seq, alt)
                        hit=(this_gene_copy.name+":"+kind)
                        if kind[0]==kind[-1]:
                            kind='S'
                        else:
                            kind='N'
                    else:
                        hit='NC'
                        kind='NC'
                fields.insert(5, hit)
                fields.insert(6, kind)
                print>>OUT, fields[0], fields[1], fields[3], fields[4], \
                    fields[5], fields[6], fields[7], fields[8]

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
    				#while (pos>min(this_gene.stop, this_gene.start) ):
    				while (pos>max(this_gene.start, this_gene.stop ) ):
    					this_gene=genes.pop(0)

    				if this_gene.strand=='+':
    					if (pos>this_gene.start):
    						kind=self.getmut(this_gene.start, pos, seq, alt)
    						hit=(this_gene.name+":"+kind)
    						if kind[0]==kind[-1]:
    							kind='S'
    						else:
    							kind='N'
    					else:
    						hit='NC'
    						kind='NC'
    				elif this_gene.strand=='-':
    					if (pos>this_gene.start):
    						kind=self.getrev(this_gene.stop, pos, seq, alt)
    						hit=(this_gene.name+":"+kind)
    						if kind[0]==kind[-1]:
    							kind='S'
    						else:
    							kind='N'
    					else:
    						hit='NC'
    						kind='NC'
    			if max(ChiPoly)>self.CUTOFF or max(ChiFixed)>self.CUTOFF:
    				print '\t'.join(map(str, fields[:6]))+'\t'+hit+'\t'+kind+'\t'+'\t'.join(map(str, Freq))
    				print>> OUT, '\t'.join(map(str, fields[:6]))+'\t'+hit+'\t'+kind+'\t'+'\t'.join(map(str, Freq))
    		else:
    			print line,
    			fields=line.split('\t')
    			if fields[0]=="@SCFNAME       ":
    				self.isopen=True
    	OUT.close

class filter:

    def __init__(self, pathIN, pathOUT):
        self.pathIN = pathIN
        self.pathOUT = pathOUT

    def getSNPTable(self):
        dfs = []
        samples = []
        for i in os.listdir(self.pathIN):
            if i.endswith(".txt"):
                sample_name = i.split('.')[0]
                sample = re.split(r'[-_]+', sample_name)[2]
                samples.append(sample)
                names = ['Scaffold', 'Pos', 'Ref', 'Alt', 'Gene', \
                    'Coding', 'Qual', 'AC']
                IN = pd.read_csv(self.pathIN + '/' + i, delimiter = ' ', header = None)
                IN.columns = names
                # remove original header and AC column
                IN = IN.iloc[1:,:-1]
                dfs.append(IN)
            else:
                continue
        merged_df = reduce(lambda left,right: pd.merge(left,right,on=['Scaffold', 'Pos','Ref','Gene'], how='outer'), dfs)
        toRename = ['Alt', 'Coding', 'Qual']
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
        merged_df.to_csv(self.pathOUT +'_GATK.txt' ,sep='\t', \
            index = False)
