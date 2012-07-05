import sys
import collections
import math
import numpy
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version: 2.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python frequency_spectrum.py -i candidates.cmh -j full.cmh -b 100 -p 1,5,9+2,6,10+3,4,8 -s 1,2 -o out -r 1
2)	This script needs the statlib package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script calculates the minor allele frequency distinbutions (-m na) or major effect allele frequency (the allele which rises in allelefreq between two timepoint:s-m 1,2 [stands for 1st and second item in -p], should be the timepoints in comparison of the cmh test) distribution of the first population in -p for two inputs (e.g. candidate dataset and full snp dataset) and calculates the ratio for a better comparison. Furthermore, it prints the average coverage for every bin in the distribution plus standard deviations. It also calculates the absolute allelefrequency changes in all populations within one file to see how much allelefrequencies are changing depending on the initial (minor) allele frequency. This results are printed in and outputfile with the ending \".distribution\". See example below:

bins	count_B_F15_Cand_1in10000	mean_cov_B_F15_Cand_1in10000	SD_cov_B_F15_Cand_1in10000	count_B_F15_Cand_1in1000	mean_cov_B_F15_Cand_1in1000	SD_cov_B_F15_Cand_1in1000	ratio
0.0	380	82.4026315789	38.1413631304	380	82.4026315789	38.1413631304	1.0
0.05	135	73.1185185185	22.3346500936	135	73.1185185185	22.3346500936	1.0
0.1	96	73.375	23.3757634269	96	73.375	23.3757634269	1.0
0.15	87	80.0459770115	22.4130921281	87	80.0459770115	22.4130921281	1.0
0.2	67	73.2089552239	19.69582121	67	73.2089552239	19.69582121	1.0
0.25	64	77.09375	24.4817513833	64	77.09375	24.4817513833	1.0
0.3	69	78.4637681159	17.321849268	69	78.4637681159	17.321849268	1.0
0.35	54	76.8148148148	25.3578718187	54	76.8148148148	25.3578718187	1.0
0.4	47	77.1276595745	23.0534511434	47	77.1276595745	23.0534511434	1.0
0.45	39	79.0	16.9876115852	39	79.0	16.9876115852	1.0


bins: 					number and width of bins
count_inp1: 			count of SNPs in corresponding bins from input1
mean_cov_inp1: 			average coverage of SNPs in input1
SD_cov_cov_inp1: 		SD of coverage of SNPs in input1
count_inp2: 			count of SNPs in corresponding bins from input2
mean_cov_inp2: 			average coverage of SNPs in input2
SD_cov_cov_inp2: 		SD of coverage of SNPs in input2
ratio_of_counts:			count of SNPs of inp1 divided by count of SNPs of inp1


4) It may happen that the number of SNPs in the bins differs from the input. If a SNP is not polymorphic in the replicate defined with -p it will be exluded.
	""") 
#########################################################   CODE   #########################################################################


parser.add_option("-i", "--inp1", dest="inp1", help="input1: cmh or sync file")
parser.add_option("-b", "--bins", dest="bins", help="No. of bins for the histogram")
parser.add_option("-p", "--pops", dest="pops", help="define, for which population in the sync you want to produce the minor allele frequency distribution")
parser.add_option("-s", "--se", dest="se", help="define start and end timepoints in analysis, e.g. Base-F15 = 1,2")
parser.add_option("-o", "--out", dest="out", help="out file")
parser.add_option("-r", "--rep", dest="rep", help="population to be tested")
parser.add_option("-c", "--cutoff", dest="cutoff", help="allele-frequency cutoff")

parser.add_option_group(group)
(options, args) = parser.parse_args()

#2L	4910	A	6:0:0:0:0:0	-	5:0:0:0:0:0	25:0:0:0:0:0	33:0:0:0:0:0
#2L	4911	G	0:0:0:7:0:0	-	0:0:0:5:0:0	0:0:0:27:0:0	0:0:0:33:0:0
#2L	4912	A	7:0:0:0:0:0	-	5:0:0:0:0:0	27:0:0:0:0:0	32:0:0:0:0:0
#2L	4913	G	0:0:0:7:0:0	-	0:0:0:5:0:0	0:0:0:27:0:0	0:0:0:34:0:0
#2L	4914	A	7:0:0:0:0:0	-	5:0:0:0:0:0	26:0:0:0:0:0	33:0:0:0:0:0
#2L	4915	G	0:0:0:8:0:0	-	0:0:0:5:0:0	0:0:0:22:0:0	0:0:0:33:0:0
#2L	4916	A	8:0:0:0:0:0	-	6:0:0:0:0:0	25:0:0:0:0:0	33:0:0:0:0:0
#2L	4917	G	0:0:0:8:0:0	-	0:0:0:6:0:0	0:0:0:23:0:0	0:0:0:32:0:0
#2L	4918	C	0:0:8:0:0:0	-	0:0:6:0:0:0	0:0:22:0:0:0	0:0:30:0:0:0
#print "###chr\tpos\trc\tallele_states\t+"""

############################################## functions ############################################################
def freqhash(inp,i):
	h,h2=collections.defaultdict(lambda:0),{}
	h["A"]+=int(inp.split()[i+2].split(":")[0])
	h["T"]+=int(inp.split()[i+2].split(":")[1])
	h["C"]+=int(inp.split()[i+2].split(":")[2])
	h["G"]+=int(inp.split()[i+2].split(":")[3])
	cov=sum(h.values())
	if cov!=0:
		for k,v in h.items():
			h2[k]=float(v)/cov
		#print h2
		return h2
	else:
		return "na"
	
def meanstdv(x): ### calculate mean, median, stdev. standard error : x=datalist
	from math import sqrt
	n,mean,std,se,median = len(x),0,0,0,0
	for a in x:
		mean = mean + a
	mean = mean / float(n)
	for a in x:
		std = std + (a - mean)**2
	std = sqrt(std / float(n-1))
	se= std/sqrt(n)
	median=sorted(x)[len(x)/2]
	return mean,std,se,median

def max_allele_change(inp,cp):
	h,h1,h2={},collections.defaultdict(lambda:0),collections.defaultdict(lambda:0)
	for i in cp[0]:
		#print i
		h1["A"]+=int(inp.split()[int(i)+2].split(":")[0])
		h1["T"]+=int(inp.split()[int(i)+2].split(":")[1])
		h1["C"]+=int(inp.split()[int(i)+2].split(":")[2])
		h1["G"]+=int(inp.split()[int(i)+2].split(":")[3])
	for i in cp[1]:
		#print i
		h2["A"]+=int(inp.split()[int(i)+2].split(":")[0])
		h2["T"]+=int(inp.split()[int(i)+2].split(":")[1])
		h2["C"]+=int(inp.split()[int(i)+2].split(":")[2])
		h2["G"]+=int(inp.split()[int(i)+2].split(":")[3])
	for k,v in h1.items():
		h[k]=h2[k]-v
	if len(set(h.values()))==1.0:
		return "na"
	else:
		maxallele=max(h,key=lambda a: h.get(a))
		return maxallele
		
def minallelefreq_mea(inp,r,p,cp):
	minl,minus=[],0
	#print p
	for l in open(inp,"r"):
		if l.rstrip()!="":
			a=l.split()
			for i in range(0,len(p[0].split(","))):
				if "-" in a[int(p[0].split(",")[i])+2]:
					minus=1
			if minus==0:
				min1=max_allele_change(l,cp)
				if min1!="na":
					minallelefr=freqhash(l,int(r))[min1]
					coverage=cov(l,int(r))
					minl.append([minallelefr,coverage,l])
				else:
					print l.rstrip()
	print len(minl)
	return minl


def allelefreq_diff(inp,cutoff,cp,out,cr):
	
	outth=open(out+"_larger-than"+str(cutoff),"w")
	for l in open(inp,"r"):
		diffl=[]
		mea=max_allele_change(l,cr)
		for i in range(0,len(cp)):
			pop1freq=freqhash(l,int(cp[i][0]))[mea]
			pop2freq=freqhash(l,int(cp[i][1]))[mea]
			diffl.append(pop2freq-pop1freq)
		if meanstdv(diffl)[0]>=cutoff:
			outth.write(l)
				
def summary_cand(inp,bin):
	ccand,covav,covsd,bins,covl=[],[],[],[],[]
	min1,c=0.0,0
	for i in range(1,int(bin)+1):
		max1=min1+(1/float(bin))
		out=open(str(options.out)+"_"+str(min1)+"-"+str(max1),"w")
		for minor,cov,l in inp:
			#print minor,i
			if minor>=min1 and minor<max1:
				c+=1
				covl.append(cov)
				out.write(l)
		bins.append(str(min1))
		ccand.append(str(c))
		if len(covl)>1:
			covav.append(str (meanstdv(covl)[0]))
			covsd.append(str(meanstdv(covl)[1]))
		if len (covl)==1:
			covav.append(covl[0])
			covsd.append(0)
		else:
			covav.append("na")
			covsd.append("na")	
			c=0
		covl=[]
		min1=max1
	return ccand,covav,covsd,bins

def cov(inp,i): ##returns coverage of population i
	h=collections.defaultdict(lambda:0)
	h["A"]=int(inp.split()[int(i)+2].split(":")[0])
	h["T"]=int(inp.split()[int(i)+2].split(":")[1])
	h["C"]=int(inp.split()[int(i)+2].split(":")[2])
	h["G"]=int(inp.split()[int(i)+2].split(":")[3])
	cov=sum(h.values())
	if cov!=0:
		return cov
	else:
		return "na"
	
############################################## code ############################################################	
rep=int(options.rep)
pops=str(options.pops).split("+")
start=int(str(options.se).split(",")[0])-1
end=int(str(options.se).split(",")[1])-1
co=float(options.cutoff)
o,allpops,comptime=[],[],[]
for line in pops:
	o.append(line.split(","))

genepops=zip(*o)
comprep=[genepops[start],genepops[end]] 																						# replicates at different timepoints
for l in comprep:	
	for i in l:
		allpops.append(i)
for l in pops:
	comptime.append((l.split(",")[start],l.split(",")[end])) 	

cand=minallelefreq_mea(str(options.inp1),rep,pops,comprep)
#candth=allelefreq_diff(str(options.inp1),co,comptime,options.out,comprep)

ccand,covav,covsd,bins=summary_cand(cand,str(options.bins))

out2=open(str(options.out)+"_mea.distributions","w")
out2.write( "bins\tcount_"+str(options.inp1).split("/")[-1:][0]+"\tmean_cov_"+str(options.inp1).split("/")[-1:][0]+"\tSD_cov_"+str(options.inp1).split("/")[-1:][0]+"\n")
new=zip(bins,ccand,covav,covsd)				
for a in new:
	out2.write( str(a[0])+"\t"+str(a[1])+"\t"+str(a[2])+"\t"+str(a[3])+"\n")
