import sys
import math
from optparse import OptionParser, OptionGroup
import collections

#Author: Martin Kapun & Robert Kofler

#########################################################   HELP   #########################################################################
usage="\npython %prog --input individuals.consensus --min-count 3 --in2lt 2,3,6,10 --in3rmo 2,3,7 --in3rc 0,1,4,6,10 --cold-hot 0,1,2,3 --in3rp 8 --all 0,1,2,3,4,6,7,8,10 --names 52,53,80,100,117,136,150,168,106 --output output_fstpi"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

The purpose of this script is to calculate FST and pi for different sets of indivduals with and without an inversion or which are belonging to a particular group (cold/hot). A consensus file is needed as an input (--input). All individuals, which are used in the analysis need to be defined with their poistion (0-base!!) in the consensus file. If you E.g. want to use indivduals 1,2,5 you need to put --all 0,1,4. The names of these indivduals need to be defined with --names. It may happen that allelic information is not available for all indivduals. Therefore, you need to define with --min-count the minimum number of individuals allowed fo calculating pi and FST. Then you need to define which individuals are carry a particular inversion. You can do this for In(2L)t, In(2R)Ns, In(3L)P, In(3R)Mo in In(3R)C. Per definition all indivduals not assigned to a particular inversion are considered not to carry this inversion. These indivduals are used to calcuate pi for non-inverted indviduals and for caluclating FST. Not that Inversions not used in the commandline are not considered. 
This script will produce one tab delimited output for FST, which contain information on the chromosme, position, inversion specific FST, used alleles and used indivduals (see header for detail).
The second output will be a tab-delimited text file containing the pi values of indivduals with and without inversions (see header)
The option chr defines if pi and FST should only be calculated for a subset of all chromosomes

""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="a consensus output file")
parser.add_option("--min-count", dest="m", help="minimum count of individuals for calculating pi for FST")
parser.add_option("--output", dest="o", help="output-file")
parser.add_option("--data1", dest="d1", help="dataset 1, e.g. inverted individuals")
parser.add_option("--data2", dest="d2", help="dataset 2, e.g. non-inverted individuals")
parser.add_option("--all", dest="a", help="positions of all indivduals used for the analysis")
parser.add_option("--names", dest="n", help="names of all indivduals used for the analysis")
parser.add_option("--chr", dest="c", help="chromosomes used for the analysis",default="all")

parser.add_option_group(group)
(options, args) = parser.parse_args()

def pi(x,n):
	''' calculate pi on a SNP-wise basis. where x is a vector of all allelefreqs and n is the samplesize)'''
	corr=(n-1)/float(n)
	freqsum=sum([y**2 for y in x])
	return (1-freqsum)*corr 

def fst(pa,pb,pt):
	if len(pa)==0 or len(pb)==0:
		return "NAN","NAN","NAN"
	pia=pi([v/float(sum(pa.values())) for (k,v) in pa.items()],sum(pa.values()))	
	pib=pi([v/float(sum(pb.values())) for (k,v) in pb.items()],sum(pb.values()))
	pit=pi([v/float(sum(pt.values())) for (k,v) in pt.items()],sum(pt.values()))
	#print pia,pa,pb,pib
	if pit==0:
		return 0.0,0.0,0.0
	Fst=(pit-((pia+pib)/2))/pit
	if Fst<0.0:
		return 0.0,pia,pib
	else:
		return Fst,pia,pib


############## define inviduals per inversion:
#--in2lt 2,3,6,10 --in3rmo 2,3,7 --in3rc 0,1,4,6,10 --cold-hot 0,1,2,3 --in3rp 8 --all 0,1,2,3,4,6,7,8,10 --names 52,53,80,100,117,136,150,168,106

chromo=options.c.split(",")
positions=map(int,options.a.split(","))
idhash=dict(zip(positions,options.n.split(",")))
data1=map(int,options.d1.split(","))
data2=map(int,options.d2.split(","))

outfst=open(options.o+".fst","w")
outpi=open(options.o+".pi","w")

## write headers
outfst.write("Chr\tpos\tFST\talleles\tindividuals\n")
outpi.write("Chr\tpos\tpi_data1\tpi_data2\n")
		
for l in open(options.input,"r"):
	if len(l.split())==1:
		continue
	chr,pos,ref,mel36,pops=l.split()[:5]

	if chromo!=["all"]:
		if chr not in chromo:
			continue
	fstlist,pilist,allelelist,poslist="","","",""
	inv=collections.defaultdict(lambda:0)
	non=collections.defaultdict(lambda:0)
	full=collections.defaultdict(lambda:0)
	inst,nost,inp,nop=[],[],[],[] 
	## loop through all indivduals:
	for item in data1:
		if pops[item]!="N":
			## count the alleles of all individuals 
			full[pops[item]]+=1
			## count the alleles of inverted individuals 
			inv[pops[item]]+=1
			## make a list of the alleles of inverted individuals 
			inst.append(pops[item])
			## make a list of the names of inverted individuals corresponding to the alleles
			inp.append(str(idhash[item]))
	for item in data2:
		if pops[item]!="N":				
			## count the alleles of all individuals 
			full[pops[item]]+=1
			## count the alleles of non-inverted individuals 
			non[pops[item]]+=1
			## make a list of the alleles of non-inverted individuals
			nost.append(pops[item])
			## make a list of the names of non-inverted individuals corresponding to the alleles
			nop.append(str(idhash[item]))

	invcov=float(sum(inv.values()))
	noncov=float(sum(non.values()))
	fullcov=float(sum(full.values()))
	# test if the number of indivduals in both groups is above the minimum-count thresholdd
	if invcov<int(options.m) or noncov<int(options.m):
		fstlist="NAN"
		pilist="NAN\tNAN"
	else:
		## calculate pi and FST and append to lists
		fstresult=fst(inv,non,full)
		fstlist=str(fstresult[0])
		pilist="\t".join(map(str,fstresult[1:]))
	allelelist="".join(inst)+"/"+"".join(nost)
	poslist=",".join(inp)+"/"+",".join(nop)
	## write to output if FST is not "NAN" or 0 in all inversions
	if fstlist!='NAN': 
		outfst.write("\t".join([chr,pos])+"\t"+fstlist+"\t"+allelelist+"\t"+poslist+"\n")
	if pilist!="NAN\tNAN": 
		outpi.write("\t".join([chr,pos])+"\t"+pilist+"\n")
	
	