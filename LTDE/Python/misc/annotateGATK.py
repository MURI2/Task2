import sys, os, copy

mydir = os.path.expanduser("~/github/Task2/LTDE")

#strains = ['KBS0703', 'KBS0705', 'KBS0706', 'KBS0710', 'KBS071', 'KBS0713', \
#    'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802', 'KBS0812']


codon_list={
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

def is_gene(line):
    try:
        if line[2]=="CDS":
            return True
        else:
            return False
    except:
        return False




def rcomp (seq):
    s=[]
    rc={'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rseq=[]
    for s in seq:
        rseq.insert(0, rc[s])
    return ''.join(rseq)

def getmut(start, bp, seq, mut):
    A=int(bp)-(int(bp)-int(start))%3
    B=int(bp)
    return codon_list[seq[A:A+3]]+str((bp-start)/3+1)+codon_list[seq[A:B]+mut+seq[B+1:A+3]]

def getrev(start, bp, seq, mut):
    bp=int(bp)+1
    start=int(start)+1

    A=int(bp)+abs(int(start)-int(bp))%3
    B=int(bp)
    return codon_list[rcomp(seq[A-3:A])]+str((start-bp)/3+1)+codon_list[rcomp(seq[A-3:B-1]+mut+seq[B:A])]

def issyn(start, bp, seq, mut):
    A=int(bp)-(int(bp)-int(start))%3
    B=int(bp)
    return codon_list[seq[A:A+3]]==codon_list[seq[A:B]+mut+seq[B+1:A+3]]

def issynrev(start, bp, seq, mut):
    bp=int(bp)+1
    start=int(start)+1

    A=int(bp)+abs(int(start)-int(bp))%3
    B=int(bp)
    return codon_list[rcomp(seq[A-3:A])]==codon_list[rcomp(seq[A-3:B-1]+mut+seq[B:A])]

def annotateGATK(species, isopen=False, CUTOFF=25.0, start=3, stop=4, name=8, strand=6):
    if species == 'KBS0812':
        GFF = open(mydir + '/data/Bacillus_test/AL009126.3.gff')
        FFN = open(mydir + '/data/Bacillus_test/AL009126.fa')
    else:
        GFF = open(mydir + '/data/2015_SoilGenomes_Annotate/' + str(species) + '/G-Chr1.gff')
        FFN = open(mydir + '/data/2015_SoilGenomes_Annotate/' + str(species) + '/G-Chr1.ffn')
    content_list = []
    for content in os.listdir(mydir + '/data/GATK/raw/' + species): # "." means current directory
        content_list.append(content)

    gene_sizes={}
    genes=[]
    for line in GFF:
        line=line.split('\t')
        if is_gene(line):
            for x in line[name].split(';'):
                if x.split('=')[0]=="gene":
                    genes.append(gene(x.split('=')[1], int(line[start]), int(line[stop]), line[strand]) )
                    gene_sizes[x.split('=')[1]]=(int(line[start])-int(line[stop]))
    seq=[]
    for line in FFN:
        if line[0]!='>':
            line=line.strip('\n')
            seq.append(line)
    seq=' '+''.join(seq)
    this_gene=genes.pop(0)

    for x in content_list:
        File = open(mydir + '/data/GATK/raw/' +   species + '/' + x)
        OUT =  open(mydir + '/data/GATK/annotate/' + str(species) + '/' + x.split('.')[0] + '_annotate.txt', 'wr')
        this_gene_copy = copy.copy(this_gene)
        genes_copy = copy.copy(genes)
        seq_copy = copy.copy(seq)
        for line in File:
            fields=line.strip().split('\t')
            if fields[0]=='CHROM':
                fields.insert( 5, 'GENE')
                fields.insert( 6, 'CODING')
                print>>OUT, fields[0], fields[1], fields[3], fields[4], \
                    fields[5], fields[6], fields[7], fields[8]
                #OUT.write( '\t'.join( fields ) )
                continue
            ref=fields[3]
            alt=fields[4]
            pos=int(fields[1])
            if len(alt)==1:
                while ((pos>max(this_gene_copy.start, this_gene_copy.stop ))  and len(genes_copy)  >0 ):
                    this_gene_copy=genes_copy.pop(0)
                if this_gene_copy.strand=='+':
                    if (pos>this_gene_copy.start):
                        kind=getmut(this_gene_copy.start, pos, seq, alt)
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
                        kind=getrev(this_gene_copy.stop, pos, seq, alt)
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
                #OUT.write( '\t'.join( fields ) )
                print>>OUT, fields[0], fields[1], fields[3], fields[4], \
                    fields[5], fields[6], fields[7], fields[8]
                #print>> OUT, '\t'.join(map(str, fields[:6]))+'\t'+hit+'\t'+kind+'\t'+'\t'.join(map(str, Freq))

        OUT.close()

#strains = ['KBS0703', 'KBS0705', 'KBS0706', 'KBS0710', 'KBS0711', 'KBS0713', \
#    'KBS0715', 'KBS0721', 'KBS0722', 'KBS0724', 'KBS0727', 'KBS0802']
strains = ['KBS0812']

for x in strains:
    annotateGATK(x)
