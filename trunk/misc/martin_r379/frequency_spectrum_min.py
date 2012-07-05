import sys
import collections
import math
import numpy
from statlib.stats import *
from statlib import pstat
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version: 2.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python frequency_spectrum_min.py -i candidates.cmh -j full.cmh -b 100 -r 1 -o out
2)	This script needs the statlib package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script calculates the minor allele frequency distinbutions (-m na)  Furthermore, it prints the average coverage for every bin in the distribution plus standard deviations. It also calculates the absolute allelefrequency changes in all populations within one file to see how much allelefrequencies are changing depending on the initial (minor) allele frequency. This results are printed in and outputfile with the ending \".distribution\". See example below:

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
parser.add_option("-j", "--inp2", dest="inp2", help="input2: cmh or sync file")
parser.add_option("-b", "--bins", dest="bins", help="No. of bins for the histogram")
parser.add_option("-r", "--rep", dest="rep", help="define which population should be analyzed")
parser.add_option("-o", "--out", dest="out", help="out file")
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
		
def minallelefreq(inp,re):
	h,h2=collections.defaultdict(lambda:0),{}
	h["A"]+=int(inp.split()[int(re)+2].split(":")[0])
	h["T"]+=int(inp.split()[int(re)+2].split(":")[1])
	h["C"]+=int(inp.split()[int(re)+2].split(":")[2])
	h["G"]+=int(inp.split()[int(re)+2].split(":")[3])
	cov=sum(h.values())
	#print l, h,cov
	if cov!=0:
		for k,v in h.items():
			if float(v)/cov!=0.0:
				h2[k]=float(v)/cov
		if min(h2.values())==1.0:
			return 0.0,cov
		else:
			return min(h2.values()),cov
	else:
		return 0.0,cov

def minallelefreq_min(inp,r):
	minl=[]
	for l in open(inp):
		a=l.split()
		if "-" not in a[int(r)+2]:
			minl.append(minallelefreq(l,r))
	return minl

def summary(inp,bin):
	ccand,covav,covsd,bins,covl=[],[],[],[],[]
	min1,c=0.0,0
	for i in range(1,int(bin)+1):
		max1=min1+(0.5/float(bin))
		for minor,cov in inp:
			#print minor,i
			if minor>=min1 and minor<max1:
				c+=1
				covl.append(cov)
		bins.append(str(min1))
		ccand.append(str(c))
		covav.append(str (lmean(covl)))
		covsd.append(str(lstdev(covl)))
		c=0
		covl=[]
		min1=max1
	return ccand,covav,covsd,bins

############################################## code ############################################################	
	
rep=int(options.rep)
	

cand=minallelefreq_min(str(options.inp1),rep)
full=minallelefreq_min(str(options.inp2),rep)

ccand,covav,covsd,bins=summary(cand,str(options.bins))
fcand,fcovav,fcovsd,bins=summary(full,str(options.bins))

out2=open(str(options.out)+"_min.distributions","w")
out2.write( "bins\tcount_"+str(options.inp1).split("/")[-1:][0]+"\tmean_cov_"+str(options.inp1).split("/")[-1:][0]+"\tSD_cov_"+str(options.inp1).split("/")[-1:][0]+"\tcount_"+str(options.inp2).split("/")[-1:][0]+"\tmean_cov_"+str(options.inp2).split("/")[-1:][0]+"\tSD_cov_"+str(options.inp2).split("/")[-1:][0]+"\tratio\n")
new=zip(bins,ccand,covav,covsd,fcand,fcovav,fcovsd)				
for a in new:
	out2.write( str(a[0])+"\t"+str(a[1])+"\t"+str(a[2])+"\t"+str(a[3])+"\t"+str(a[4])+"\t"+str(a[5])+"\t"+str(a[6])+"\t"+str(float(a[1])/float(a[4]))+"\n")
