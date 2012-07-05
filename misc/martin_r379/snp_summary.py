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
parser.add_option("-v", "--overlap", dest="o", help="prints overlap of SNPs and genes from both inputs; put \"yes\", or \"no\"")
parser.add_option("-g", "--genelist1", dest="g", help="prints genelist for input1; put \"yes\", or \"no\"")
parser.add_option("-f", "--genelist2", dest="f", help="prints genelist for input2; put \"yes\", or \"no\"")
parser.add_option("-o", "--output", dest="p", help="define general output file")
parser.add_option_group(group)
(options, args) = parser.parse_args()

### inversion 3R(Payne) based on 
#proximal breakpoint: somewhere around 12253902
#distal breakpoint: somewhere around 20565305
fdict,cdict={},{}
def summary(data):
	genes=open(data+".genes","r")
	igv=open(data+".igv","r")
	payne_prox=int(12253902)
	payne_dist=int(20565305)
	chromosomes=["2L","2LHet","2R","2RHet","3L","3LHet","3R","3Rinv","3Rout","3RHet","4","X","Xhet","Yhet","Aut","total"]
	features=["SYNONYMOUS_STOP","STOP_LOST","NON_SYNONYMOUS_CODING","START_GAINED","3PRIME_UTR","SPLICE_SITE","INTRON","START_LOST","UPSTREAM","INTERGENIC","DOWNSTREAM","SYNONYMOUS_CODING","5PRIME_UTR","STOP_GAINED","total"]
	genelist,posdict={},{}
	fdict=collections.defaultdict(lambda:0)
	cdict=collections.defaultdict(lambda:0)
	for a in features:
		fdict[a]=0
	for a in chromosomes:
		cdict[a]=0
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
	

fdict_full,cdict_full,genelist_full,pos_full=summary(str(options.i))
fdict_cand,cdict_cand,genelist_cand,pos_cand=summary(str(options.j))


for k,v in fdict_cand.items():
	fdict[k]=[v,fdict_full[k]]
for k,v in cdict_cand.items():
	cdict[k]=[v,cdict_full[k]]

outg=open(str(options.p)+".summary","w")

if "/" in str(options.i) and "/" in str(options.j):
	outg.write("Snp_Effect\tcount:"+str(options.j).split("/")[-1]+"\tfreq:"+str(options.j).split("/")[-1]+"\tcount:"+str(options.i).split("/")[-1]+"\tfreq:"+str(options.i).split("/")[-1]+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	for k,v in sorted(fdict.items()):
		if k!="total":
			outg.write(k+"\t"+str(v[0])+"\t"+str(float(v[0])/float(fdict["total"][0]))+"\t"+str(v[1])+"\t"+str(float(v[1])/float(fdict["total"][1]))+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("total\t"+str(fdict["total"][0])+"\t"+str(float(fdict["total"][0])/float(fdict["total"][0]))+"\t"+str(fdict["total"][1])+"\t"+str(float(fdict["total"][1])/float(fdict["total"][1]))+"\n")
	outg.write("\n")
	outg.write("Chromosome\tcount:"+str(options.j).split("/")[-1]+"\tfreq:"+str(options.j).split("/")[-1]+"\tcount:"+str(options.i).split("/")[-1]+"\tfreq:"+str(options.i).split("/")[-1]+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	for k,v in sorted(cdict.items()):
		if k!="total" and k!="3Rinv" and k!="3Rout"and k!="Aut":
			outg.write(k+"\t"+str(v[0])+"\t"+str(float(v[0])/float(cdict["total"][0]))+"\t"+str(v[1])+"\t"+str(float(v[1])/float(cdict["total"][1]))+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	
	outg.write("total\t"+str(cdict["total"][0])+"\t"+str(float(cdict["total"][0])/float(cdict["total"][0]))+"\t"+str(cdict["total"][1])+"\t"+str(float(cdict["total"][1])/float(cdict["total"][1]))+"\n")
	outg.write("\n")
	outg.write("""SNPs inside or outside of In(3R)Payne:"""+"\n")
	outg.write("location\tcount:"+str(options.j).split("/")[-1]+"\tfreq:"+str(options.j).split("/")[-1]+"\tcount:"+str(options.i).split("/")[-1]+"\tfreq:"+str(options.i).split("/")[-1]+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("3Rinv\t"+str(cdict["3Rinv"][0])+"\t"+str(float(cdict["3Rinv"][0])/float(cdict["3R"][0]))+"\t"+str(cdict["3Rinv"][1])+"\t"+str(float(cdict["3Rinv"][1])/float(cdict["3R"][1]))+"\n")
	outg.write("3Rout\t"+str(cdict["3Rout"][0])+"\t"+str(float(cdict["3Rout"][0])/float(cdict["3R"][0]))+"\t"+str(cdict["3Rout"][1])+"\t"+str(float(cdict["3Rout"][1])/float(cdict["3R"][1]))+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("3R\t"+str(cdict["3R"][0])+"\t"+str(float(cdict["3R"][0])/float(cdict["3R"][0]))+"\t"+str(cdict["3R"][1])+"\t"+str(float(cdict["3R"][1])/float(cdict["3R"][1]))+"\n")
	outg.write("\n")
	outg.write("""Ratio X to autosome:"""+"\n")
	outg.write("location\tcount:"+str(options.j).split("/")[-1]+"\tfreq:"+str(options.j).split("/")[-1]+"\tcount:"+str(options.i).split("/")[-1]+"\tfreq:"+str(options.i).split("/")[-1]+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("X\t"+str(cdict["X"][0])+"\t"+str(float(cdict["X"][0])/float(cdict["total"][0]))+"\t"+str(cdict["X"][1])+"\t"+str(float(cdict["X"][1])/float(cdict["total"][1]))+"\n")
	outg.write("Aut\t"+str(cdict["Aut"][0])+"\t"+str(float(cdict["Aut"][0])/float(cdict["total"][0]))+"\t"+str(cdict["Aut"][1])+"\t"+str(float(cdict["Aut"][1])/float(cdict["total"][1]))+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("total\t"+str(cdict["total"][0])+"\t"+str(float(cdict["total"][0])/float(cdict["total"][0]))+"\t"+str(cdict["total"][1])+"\t"+str(float(cdict["total"][1])/float(cdict["total"][1]))+"\n")
else:
	outg.write("Snp_Effect\tcount:"+str(options.j)+"\tfreq:"+str(options.j)+"\tcount:"+str(options.i)+"\tfreq:"+str(options.i)+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	for k,v in sorted(fdict.items()):
		if k!="total":
			outg.write(k+"\t"+str(v[0])+"\t"+str(float(v[0])/float(fdict["total"][0]))+"\t"+str(v[1])+"\t"+str(float(v[1])/float(fdict["total"][1]))+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("total\t"+str(fdict["total"][0])+"\t"+str(float(fdict["total"][0])/float(fdict["total"][0]))+"\t"+str(fdict["total"][1])+"\t"+str(float(fdict["total"][1])/float(fdict["total"][1]))+"\n")
	outg.write("\n")
	outg.write("Chromosome\tcount:"+str(options.j)+"\tfreq:"+str(options.j)+"\tcount:"+str(options.i)+"\tfreq:"+str(options.i)+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	for k,v in sorted(cdict.items()):
		if k!="total" and k!="3Rinv" and k!="3Rout"and k!="Aut":
			outg.write(k+"\t"+str(v[0])+"\t"+str(float(v[0])/float(cdict["total"][0]))+"\t"+str(v[1])+"\t"+str(float(v[1])/float(cdict["total"][1]))+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	
	outg.write("total\t"+str(cdict["total"][0])+"\t"+str(float(cdict["total"][0])/float(cdict["total"][0]))+"\t"+str(cdict["total"][1])+"\t"+str(float(cdict["total"][1])/float(cdict["total"][1]))+"\n")
	outg.write("\n")
	outg.write("""SNPs inside or outside of In(3R)Payne:"""+"\n")
	outg.write("location\tcount:"+str(options.j)+"\tfreq:"+str(options.j)+"\tcount:"+str(options.i)+"\tfreq:"+str(options.i)+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("3Rinv\t"+str(cdict["3Rinv"][0])+"\t"+str(float(cdict["3Rinv"][0])/float(cdict["3R"][0]))+"\t"+str(cdict["3Rinv"][1])+"\t"+str(float(cdict["3Rinv"][1])/float(cdict["3R"][1]))+"\n")
	outg.write("3Rout\t"+str(cdict["3Rout"][0])+"\t"+str(float(cdict["3Rout"][0])/float(cdict["3R"][0]))+"\t"+str(cdict["3Rout"][1])+"\t"+str(float(cdict["3Rout"][1])/float(cdict["3R"][1]))+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("3R\t"+str(cdict["3R"][0])+"\t"+str(float(cdict["3R"][0])/float(cdict["3R"][0]))+"\t"+str(cdict["3R"][1])+"\t"+str(float(cdict["3R"][1])/float(cdict["3R"][1]))+"\n")
	outg.write("\n")
	outg.write("""Ratio X to autosome:"""+"\n")
	outg.write("location\tcount:"+str(options.j)+"\tfreq:"+str(options.j)+"\tcount:"+str(options.i)+"\tfreq:"+str(options.i)+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("X\t"+str(cdict["X"][0])+"\t"+str(float(cdict["X"][0])/float(cdict["total"][0]))+"\t"+str(cdict["X"][1])+"\t"+str(float(cdict["X"][1])/float(cdict["total"][1]))+"\n")
	outg.write("Aut\t"+str(cdict["Aut"][0])+"\t"+str(float(cdict["Aut"][0])/float(cdict["total"][0]))+"\t"+str(cdict["Aut"][1])+"\t"+str(float(cdict["Aut"][1])/float(cdict["total"][1]))+"\n")
	outg.write("__________________________________________________________________________________________"+"\n")
	outg.write("total\t"+str(cdict["total"][0])+"\t"+str(float(cdict["total"][0])/float(cdict["total"][0]))+"\t"+str(cdict["total"][1])+"\t"+str(float(cdict["total"][1])/float(cdict["total"][1]))+"\n")
if str(options.o)=="yes":
	outo=open(str(options.p)+".overlap","w")
	outos=open(str(options.p)+".snpoverlap","w")
	c=0
	for k,v in sorted(genelist_full.items()):
		if k in genelist_cand:
			c+=1
	outo.write( str(c)+" genes overlap between "+str(options.j)+" and "+str(options.i)+"\n")
	for k,v in sorted(genelist_full.items()):
		if k in genelist_cand:
			outo.write( k+"\t"+v+"\n")
	for l in open(str(options.j),"r"):
		if l.split()[0]+"_"+l.split()[2] in pos_full:
			outos.write(l)
if str(options.g)=="yes":
	outgl1=open(str(options.p)+".gl1","w")
	outgl1.write( "SNPs overlap with "+str(len(genelist_full))+" genes in "+str(options.i)+"\n")
	for k,v in sorted(genelist_full.items()):
		outgl1.write( k+"\t"+v+"\n")
if str(options.f)=="yes":
	outgl2=open(str(options.p)+".gl2","w")
	outgl2.write( "SNPs overlap with "+str(len(genelist_cand))+" genes in "+str(options.j)+"\n")
	for k,v in sorted(genelist_cand.items()):
		outgl2.write( k+"\t"+v+"\n")