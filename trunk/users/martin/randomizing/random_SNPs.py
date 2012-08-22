import sys
import collections
from optparse import OptionParser,OptionGroup
import random

#Author: Martin Kapun

#########################################################   HELP   #########################################################################

usage= "python %prog --inp input.sync --cand candidates.sync --samples 10 --number 1000  --output random_cand"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,"""
H E L P:
____________

This script samples n times (--samples) a number (--number) of SNPs  without replacement from the input (--input) an writes separate files for each sample. Optionally it is possible to define a set of SNPs to be excluded (--cand)
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-i", "--inp", dest="inp", help="input")
parser.add_option("-c", "--cand", dest="cand", help="candidate data to be excluded",default="no")
parser.add_option("-s", "--samples", dest="samples", help="number of samples")
parser.add_option("-n", "--number", dest="number", help="number of SNPs per sample")
parser.add_option("-o", "--out", dest="out", help="outputfile")
parser.add_option_group(group)
(options, args) = parser.parse_args()


inputhash={}
candhash={}

if options.cand!="no":
	for l in open(str(options.cand),"r"):
		candhash[l]=0

for l in open(str(options.inp),"r"):
	if l not in candhash:
		inputhash[l]=0

for i in range(int(options.samples)):
	if len(inputhash)>=int(options.number):
		out=open(str(options.out),"w")
		subsample=random.sample(inputhash,int(options.number))

		for k in subsample:
			del inputhash[k]
			out.write(k)
	else: 
		print "dataset too small for subsampling, only "+str(i)+" sampling(s) possible and already performed"
		break
		