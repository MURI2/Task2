import annotateFuncs as af
import os

mydir = os.path.expanduser("~/github/Task2/LTDE")

strains = ['KBS0812']


#for strain in strains:
#    content_list = []
#    for content in os.listdir(mydir + '/data/GATK/raw/' + strain): # "." means current directory
#        content_list.append(content)
#    GFF = mydir + '/data/Bacillus_test/AL009126.3.gff'
#    FFN = mydir + '/data/Bacillus_test/AL009126.fa'
#    mapgd=mydir + '/data/mapgd/raw/' + str(strain) + '_merged.pol'
#    mapgsOUT =  mydir+'/data/mapgd/annotate/' + str(strain) +  '_merged_annotate.pol'
#    af.annotate(GFF, FFN).annotateMAPGD(mapgd, mapgsOUT)
#    for x in content_list:
#        IN = mydir + '/data/GATK/raw/' +  strain + '/' + x
#        OUT =  mydir + '/data/GATK/annotate/' + str(strain) + '/' + x.split('.')[0] + '_annotate.txt'
#        af.annotate(GFF, FFN).annotateGATK(IN, OUT)

for strain in strains:
    pathin = mydir + '/data/GATK/annotate/' + strain
    pathout = mydir + '/data/GATK/merged/all/'  + strain
    af.filter(pathin, pathout).getSNPTable()
