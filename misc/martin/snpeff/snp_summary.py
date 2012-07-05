import sys
from optparse import OptionParser,OptionGroup
import collections 

#version 1.0
#Author: Martin Kapun

#########################################################   HELP   #########################################################################
print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python snp_summary.py -i SNPs -j candidates -o out -v yes -g yes -f yes 
2)	This script compares two snpeff outputs and produces several summaries. The script requires files with the exxtensions \".genes\" and \".igv\" for each dataset to be located in the same folder with the same name. Therefore, for parameters -i and -j put the candidate files without extensions, the program should find them . The   In the general output (out.summary), it produces tables showing the SNP counts and frequencies for (I) SNP effects and (II) distribution on the chromosomes for both inputs (-i and -j) and (III) the proportion of SNPs on the 3R chromosome located inside or outside the Payne inversion and (IV) the ratio of SNPs on the X and and the autosomes. Optionally it produces tables of genes in file 1 (-g) or file 2 (-f). Additionally the script produces a list of genes which are overlapping between the two input files (-v) and another output file containing the SNPs overlapping in both datasets.
	""")

#########################################################   CODE   #########################################################################

parser.add_option("-i", "--input1", dest="i", help="input file 1: Needs to be an output from the candidates_link_snpeff.py program")
parser.add_option("-j", "--input2", dest="j", help="input file 2: Needs to be an output from the candidates_link_snpeff.py program")
parser.add_option("-g", "--genelist2", dest="g", help="prints genelist for input2; put \"yes\", or \"no\"")
parser.add_option("-f", "--genelist1", dest="f", help="prints genelist for input1; put \"yes\", or \"no\"")
parser.add_option("-o", "--output", dest="p", help="define general output file")
parser.add_option_group(group)
(options, args) = parser.parse_args()

### inversion 3R(Payne) based on 
#proximal breakpoint: somewhere around 12253902
#distal breakpoint: somewhere around 20565305
def dict(inp1,inp2):
	dict=collections.defaultdict(lambda:0)
	for l in open(inp1+".genes","r"):
		if "Chromosome" not in l.split()[0] and "na" not in l.split()[0] and l.rstrip()!="": 
			a=l.split("\t")
			dict[a[-4]]
	for l in open(inp2+".genes","r"):
		if "Chromosome" not in l.split()[0] and "na" not in l.split()[0] and l.rstrip()!="": 
			a=l.split("\t")
			dict[a[-4]]
	return dict				
		


def summary(data,fdict):
	genes=open(data+".genes","r")
	igv=open(data+".igv","r")
	payne_prox=int(12253902)
	payne_dist=int(20565305)
	chromosomes=["2L","2LHet","2R","2RHet","3L","3LHet","3R","3Rinv","3Rout","3RHet","4","X","Xhet","Yhet","Aut","total"]
	genelist,posdict={},{}
	cdict={}.fromkeys(chromosomes,0)
	for l in igv:
		if "Chromosome" not in l.split()[0] and "na" not in l.split()[0] and l.rstrip()!="":
			chrom=l.split("\t")[0]
			pos=l.split("\t")[2]
			cdict["total"]+=1
			cdict[chrom]+=1
			if chrom=="3R" and int(pos)>payne_prox and int(pos)<payne_dist:
				cdict["3Rinv"]+=1
			if chrom=="3R" and int(pos)<payne_prox:
				cdict["3Rout"]+=1
			if chrom=="3R" and int(pos)>payne_dist:
				cdict["3Rout"]+=1
			if chrom!="X":
				cdict["Aut"]+=1
	for l in genes:
		if "Chromosome" not in l.split()[0] and "na" not in l.split()[0] and l.rstrip()!="":
			posdict[l.split()[0]+"_"+l.split()[2]]=1
			fdict["total"]+=1
			chrom=l.split("\t")[0]
			pos=l.split("\t")[2]
			effect=l.split("\t")[-4:-3]
			gene_name=l.split("\t")[-5:-4]
			gene_id=l.split("\t")[-6:-5]
			if "na" not in gene_name:
				genelist[gene_name[0]]=gene_id[0]
			fdict[effect[0]]+=1
			
			
	return fdict,cdict,genelist,posdict

def names(inp):
	if "/" in inp:
		name=inp.split("/")[-1]
	else:
		name=inp
	return name

#print dict(str(options.i),str(options.j))
	
fdict_full,cdict_full,genelist_full,pos_full=summary(str(options.j),dict(str(options.i),str(options.j)))
fdict_cand,cdict_cand,genelist_cand,pos_cand=summary(str(options.i),dict(str(options.i),str(options.j)))

outg=open(str(options.p)+".summary","w")

names_i=names(str(options.i))
names_j=names(str(options.j))
#### print summary of SNP effects ######
outg.write("Snp_Effect\tcount:"+names_j+"\tfreq:"+names_j+"\tcount:"+names_i+"\tfreq:"+names_i+"\n")
outg.write("__________________________________________________________________________________________"+"\n")
for effect,v in sorted(fdict_full.items()):
	if effect!="total":
		outg.write(effect+"\t"+str(v)+"\t"+str(float(v)/float(fdict_full["total"]))+"\t"+str(fdict_cand[effect])+"\t"+str(float(fdict_cand[effect])/float(fdict_cand["total"]))+"\n")
outg.write("__________________________________________________________________________________________"+"\n")
outg.write("total\t"+str(fdict_full["total"])+"\t"+str(float(fdict_full["total"])/float(fdict_full["total"]))+"\t"+str(fdict_cand["total"])+"\t"+str(float(fdict_cand["total"])/float(fdict_cand["total"]))+"\n")
outg.write("\n")
#### print summary of distribution across genome ######
outg.write("Chromosome\tcount:"+names_j+"\tfreq:"+names_j+"\tcount:"+names_i+"\tfreq:"+names_i+"\n")
outg.write("__________________________________________________________________________________________"+"\n")
for chromosome,v in sorted(cdict_full.items()):
	if chromosome!="total" and chromosome!="3Rinv" and chromosome!="3Rout"and chromosome!="Aut":
		outg.write(chromosome+"\t"+str(v)+"\t"+str(float(v)/float(cdict_full["total"]))+"\t"+str(cdict_cand[chromosome])+"\t"+str(float(cdict_cand[chromosome])/float(cdict_cand["total"]))+"\n")
outg.write("__________________________________________________________________________________________"+"\n")
outg.write("total\t"+str(cdict_full["total"])+"\t"+str(cdict_full["total"]/float(cdict_full["total"]))+"\t"+str(cdict_cand["total"])+"\t"+str(float(cdict_cand["total"])/float(cdict_cand["total"]))+"\n")
outg.write("\n")
outg.write("""SNPs_inside_or_outside_of_In(3R)Payne:"""+"\n")
outg.write("location\tcount:"+names_j+"\tfreq:"+names_j+"\tcount:"+names_i+"\tfreq:"+names_i+"\n")
outg.write("__________________________________________________________________________________________"+"\n")
outg.write("3Rinv\t"+str(cdict_full["3Rinv"])+"\t"+str(float(cdict_full["3Rinv"])/float(cdict_full["3R"]))+"\t"+str(cdict_cand["3Rinv"])+"\t"+str(float(cdict_cand["3Rinv"])/float(cdict_cand["3R"]))+"\n")
outg.write("3Rout\t"+str(cdict_full["3Rout"])+"\t"+str(float(cdict_full["3Rout"])/float(cdict_full["3R"]))+"\t"+str(cdict_cand["3Rout"])+"\t"+str(float(cdict_cand["3Rout"])/float(cdict_cand["3R"]))+"\n")
outg.write("__________________________________________________________________________________________"+"\n")
outg.write("3R\t"+str(cdict_full["3R"])+"\t"+str(float(cdict_full["3R"])/float(cdict_full["3R"]))+"\t"+str(cdict_cand["3R"])+"\t"+str(float(cdict_cand["3R"])/float(cdict_cand["3R"]))+"\n")
outg.write("\n")
outg.write("""Ratio X to autosome:"""+"\n")
outg.write("location\tcount:"+names_j+"\tfreq:"+names_j+"\tcount:"+names_i+"\tfreq:"+names_i+"\n")
outg.write("__________________________________________________________________________________________"+"\n")
outg.write("X\t"+str(cdict_full["X"])+"\t"+str(float(cdict_full["X"])/float(cdict_full["total"]))+"\t"+str(cdict_cand["X"])+"\t"+str(float(cdict_cand["X"])/float(cdict_cand["total"]))+"\n")
outg.write("Aut\t"+str(cdict_full["Aut"])+"\t"+str(float(cdict_full["Aut"])/float(cdict_full["total"]))+"\t"+str(cdict_cand["Aut"])+"\t"+str(float(cdict_cand["Aut"])/float(cdict_cand["total"]))+"\n")
outg.write("__________________________________________________________________________________________"+"\n")
outg.write("total\t"+str(cdict_full["total"])+"\t"+str(float(cdict_full["total"])/float(cdict_full["total"]))+"\t"+str(cdict_cand["total"])+"\t"+str(float(cdict_cand["total"])/float(cdict_cand["total"]))+"\n")

if str(options.g)=="yes":
	outgl1=open(str(options.p)+".gl2","w")
	outgl1.write( "SNPs overlap with "+str(len(genelist_full))+" genes in "+str(options.j)+"\n")
	for k,v in sorted(genelist_full.items()):
		outgl1.write( k+"\t"+v+"\n")
if str(options.f)=="yes":
	outgl2=open(str(options.p)+".gl1","w")
	outgl2.write( "SNPs overlap with "+str(len(genelist_cand))+" genes in "+str(options.i)+"\n")
	for k,v in sorted(genelist_cand.items()):
		outgl2.write( k+"\t"+v+"\n")