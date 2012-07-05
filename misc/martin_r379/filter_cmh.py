import sys
import collections
from optparse import OptionParser, OptionGroup
from optparse import OptionParser,OptionGroup
import copy

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python filter_cmh.py -y candidates_Base-F15.cmh -p 1,5,9,2,6,10,3,4,8 -t 2
2)	This script filters the final CMH for consistency in alleles which change in frequency across all replicates in the comparisons (e.g. \"A\" in all three replicates). and tests whether there are differences in the polymorphism of SNPs across replicates and removes a SNP if the SNP is monomorphic in one of the selected populations. -y defines the cmh-test output, which should be filtered, -p tells the program the structure of the timepoints and replicates in the dataset: a \"+\" separates reploicates and a \",\" separates timepoints in a timeline (e.g. 1,2+3,4 would mean: 2 replicates (1,2 and 3,4) with two time points (1 and 2 for replicate 1 and 3 and 4 for replicate 2) -s defines which timepoints were used for the comparison in the CMH test of the input file: e.g if you have three timepoints (Base, F15,F37) and compared timepoint 1 with timepoint 3 (Base vs. F37), you have to put \"-s 1,3\". With the parameter -f you can further filter for consistency in MEA (the allele rising the most between two timepoints (Major effect allele: MEA)) between two replicates: e.g. -f 1,2 would compare the MEA between replicate 1 and 2 and discard the SNP if the NMEA is different. You can disable this parameter if you put \"na\" (e.g -f 1,3)
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-y", "--sync", dest="sync", help="*.sync or cmh output file")
parser.add_option("-p", "--pops", dest="pops", help="define replicates and populations: separate replicates with a \"+\" and populations with a \",\", e.g. 1,5,9,2,6,10,3,4,8 ")
parser.add_option("-t", "--th", dest="th", help="min count threshold")
parser.add_option("-o", "--out", dest="out", help="output file")

parser.add_option_group(group)
(options, args) = parser.parse_args()


def freqhash(inp,i):
	h,h2=collections.defaultdict(lambda:0),{}
	h["A"]+=int(inp.split()[i+2].split(":")[0])
	h["T"]+=int(inp.split()[i+2].split(":")[1])
	h["C"]+=int(inp.split()[i+2].split(":")[2])
	h["G"]+=int(inp.split()[i+2].split(":")[3])
	cov=sum(h.values())
	if cov!=0:
		for k,v in h.items():
			if v>=int(options.th):
				h2[k]=float(v)/cov
		if len(h2)>1:
			return 1
		else:
			return 0
	else:
		return 0

pop=str(options.pops).split(",")
out=open(str(options.out)+"_discarded","w")
out1=open(str(options.out)+"_filtered","w")
snpc=0
for l in open(str(options.sync),"r"):
	snpc+=1
	if l.rstrip()!="":
		a=l.split()
		minus=0
		for i in pop:
			if "-" in a[int(i)+2]:
				minus=1
		if minus==0:
			count=0
			for i in pop:
				count+=freqhash(l,int(i))
			if count!=1:
				out1.write(l)
			else: 
				out.write(l)
	if snpc%10000==0:
		print str(snpc)+" SNPs processed"
