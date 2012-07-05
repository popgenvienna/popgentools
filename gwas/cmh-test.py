#!/usr/bin/env python
import sys
import collections
import modules.RCMH
from modules.CMH import SyncReader
from optparse import OptionParser, OptionGroup


#Author: Robert Kofler


#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,
"""
    			D E S C R I P T I O N
Performs the CMH-test
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--output", dest="output", help="A output file")
parser.add_option("--min-count",dest="mincount",help="The minimum allele count")
parser.add_option("--min-coverage",dest="mincoverage",help="The minimum coverage")
parser.add_option("--max-coverage",dest="maxcoverage",help="The maximum coverage")
parser.add_option("--population",dest="populations",help="Calculate the p-value for the given populations")
parser.add_option("--min-pvalue",dest="minpvalue",help="The minimum p-value")
parser.add_option("--test",action="store_true", dest="test",help="run the doctest")
parser.add_option_group(group)
(options, args) = parser.parse_args()


mincount=int(options.mincount)
mincoverage=int(options.mincoverage)
maxcoverage=int(options.maxcoverage)
minpvalue=float(options.minpvalue)
popstotest=map(int,options.populations.split(","))
if not mincount:
	raise Exception("No min value")


ofw=open(options.output,"w")

for sync in SyncReader(options.input):
	
	# extract the correct populations
	pops=sync.subpopulations(popstotest)
	# test if coverage is valid
	if not modules.RCMH.Utility.coverageValid(pops, mincoverage, maxcoverage):
		continue
	
	# exclude deletions
	delcount=modules.RCMH.Utility.getDeletionCount(pops)
	if delcount>=1:
		continue
	
	# obtain counts for the two major alleles
	(alcount,majora,minora)=modules.RCMH.Utility.getMajorAlleleCount(pops)
	
	# check if SNP
	if not modules.RCMH.Utility.issnp(alcount,mincount):
		continue
	
	#perform cochran maentel haenszel test
	cmhpvalue = modules.RCMH.Utility.getCMHPvalue(alcount)
	
	# discard values which are to high
	if cmhpvalue> minpvalue:
		continue
	
	sync.pvalue=cmhpvalue
	ofw.write(str(sync)+"\n")
	

	
    