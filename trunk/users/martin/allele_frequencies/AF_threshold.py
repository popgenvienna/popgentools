import sys
import collections
import math
import numpy
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input candidates_BE.cmh --pops 8,9,10 --threshold 0.95 > candidates_BE_095.cmh"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,"""	

H E L P :
_________

Description:
This script extracts the SNP with a major allele equal or higher than the threshold (--threshold) in all replicate populations (defined with --pops) from a sync input file (--input)
""")
#########################################################   CODE   #########################################################################


parser.add_option("-i", "--input", dest="inp", help="input: cmh or sync file")
parser.add_option("-p", "--pops", dest="pops", help="define, for which population in the sync you want to test for fixed SNPs")
parser.add_option("-t", "--threshold", dest="th", help="fixed SNP threshold: e.g: 0.9")

parser.add_option_group(group)
(options, args) = parser.parse_args()

############################################## functions ############################################################
def freqhash(inp,i):
	''' return the allelefrequencies in a dictionary'''
	h,h2=collections.defaultdict(lambda:0),{}
	h["A"]+=int(inp.split()[i+2].split(":")[0])
	h["T"]+=int(inp.split()[i+2].split(":")[1])
	h["C"]+=int(inp.split()[i+2].split(":")[2])
	h["G"]+=int(inp.split()[i+2].split(":")[3])
	cov=sum(h.values())
	if cov!=0:
		for k,v in h.items():
			h2[k]=float(v)/cov
		return h2
	else:
		return "na"

	
for l in open(options.inp,"r"):
	a=l.split()
	flag=0
	for i in options.pops.split(","):
		if max(freqhash(l,int(i)).values())>float(options.th):
		        flag+=1
	if flag==len(options.pops.split(",")):
		print l.rstrip()