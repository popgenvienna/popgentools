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
The second output will be a tab-delimited text file containing the pi values of indivduals with and without inversions (see headerI
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="a consensus output file")
parser.add_option("--min-count", dest="m", help="minimum count of individuals for calculating pi for FST")
parser.add_option("--output", dest="o", help="output-file")
parser.add_option("--in2lt", dest="In2Lt", help="position of individuals carrying In(2L)t, 0-based and comma separated, per default: 'NA' ",default="NA")
parser.add_option("--in3lp", dest="In3LP", help="position of individuals carrying In(3L)P, 0-based and comma separated, per default: 'NA' ",default="NA")
parser.add_option("--in2rns", dest="In2RNs", help="position of individuals carrying In(2R)Ns 0-based and comma separated, per default: 'NA' ",default="NA")
parser.add_option("--in3rp", dest="In3RP", help="position of individuals carrying In(3R)P 0-based and comma separated, per default: 'NA' ",default="NA")
parser.add_option("--in3rmo", dest="In3RMo", help="position of individuals carrying In(3R)Mo 0-based and comma separated, per default: 'NA' ",default="NA")
parser.add_option("--in3rc", dest="In3RC", help="position of individuals carrying In(3R)C 0-based and comma separated, per default: 'NA' ",default="NA")
parser.add_option("--cold-hot", dest="coldhot", help="position of individuals carrying In(3R)C 0-based and comma separated, per default: 'NA' ",default="NA")
parser.add_option("--all", dest="a", help="positions of all indivduals used for the analysis")
parser.add_option("--names", dest="n", help="names of all indivduals used for the analysis")

parser.add_option_group(group)
(options, args) = parser.parse_args()

def pi(x,n):
	''' calculate pi on a SNP-wise basis. where x is a vector of all allelefreqs and n is the samplesize)'''
	if x=="nan":
		return "nan"
	else:
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

inversions=collections.defaultdict(lambda:"")
groups=[options.In2Lt,options.In2RNs,options.In3LP,options.In3RMo,options.In3RP,options.In3RC,options.coldhot]
types=["In2Lt","In2RNs","In3LP","In3RMo","In3RP","In3RC","cold/hot"]
positions=map(int,options.a.split(","))
idhash=dict(zip(positions,options.n.split(",")))


for i in range(len(groups)):
	if groups[i]=="NA":
		continue
	inversions[types[i]]=map(int,groups[i].split(","))


### Payne is a special case in my dataset, therefore I need to exclude it when analysing the other inversions on 3R
if "In3RP" in inversions:
	payne=inversions["In3RP"][0]	
else:
	payne="NA"


outfst=open(options.o+".fst","w")
outpi=open(options.o+".pi","w")

## write headers
outfst.write("Chr\tpos\t"+"\t".join(sorted(inversions.keys()))+"\t"+"\t".join(["alleles:"+x for x in sorted(inversions.keys())])+"\t"+"\t".join(["individuals:"+x for x in sorted(inversions.keys())])+"\n")
outpi.write("Chr\tpos\t"+"\t".join(["pi_inv:"+x+"\tpi_out:"+x for x in sorted(inversions.keys())])+"\n")
		
for l in open(options.input,"r"):
	chr,pos,ref,mel36,pops,qual,alleles=l.split()
	fstlist,pilist,allelelist,poslist=[],[],[],[]
	## loop through all inversions
	for inversion,pop in sorted(inversions.items()):
		inv=collections.defaultdict(lambda:0)
		non=collections.defaultdict(lambda:0)
		full=collections.defaultdict(lambda:0)
		inst,nost,inp,nop=[],[],[],[] 
		## loop through all indivduals:
		for item in range(len(pops)):
			## exclude In(3R)Payne haploytpe from the other 3R inversions
			if item==payne: 
				if inversion=="In3RC" or inversion=="In3RMo":
					continue
			
			if pops[item]!="N":
			
				if item in pop:
					## count the alleles of all individuals 
					full[pops[item]]+=1
					## count the alleles of inverted individuals 
					inv[pops[item]]+=1
					## make a list of the alleles of inverted individuals 
					inst.append(pops[item])
					## make a list of the names of inverted individuals corresponding to the alleles
					inp.append(str(idhash[item]))
				else:
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
			fstlist.append("NAN")
			pilist.append("NAN\tNAN")
		else:
			## calculate pi and FST and append to lists
			fstresult=fst(inv,non,full)
			fstlist.append(str(fstresult[0]))
			pilist.append("\t".join(map(str,fstresult[1:])))
		allelelist.append("".join(inst)+"/"+"".join(nost))
		poslist.append(",".join(inp)+"/"+",".join(nop))
	## write to output if FST is not "NAN" or 0 in all inversions
	if list(set(fstlist))!=['NAN'] and list(set(fstlist))!=['0.0','NAN'] and list(set(fstlist))!=['0.0']:	
		#print list(set(fstlist))	
		outfst.write("\t".join([chr,pos])+"\t"+"\t".join(fstlist)+"\t"+"\t".join(allelelist)+"\t"+"\t".join(poslist)+"\n")
		outpi.write("\t".join([chr,pos])+"\t"+"\t".join(pilist)+"\n")
	
	