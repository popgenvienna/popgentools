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
Filters SNP positions from a synchronized file
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--population",dest="populations",help="Determine whether the given position is a SNP based on these populations")
parser.add_option("--test",action="store_true", dest="test",help="run the doctest")
parser.add_option_group(group)
(options, args) = parser.parse_args()

popstotest=map(int,options.populations.split(","))

for sync in SyncReader(options.input):
	
	# extract the correct populations
	pops=sync.subpopulations(popstotest)
	
	# obtain counts for the two major alleles
	(alcount,majora,minora)=modules.RCMH.Utility.getMajorAlleleCount(pops)
	
	toprint=[]
	toprint.append("data")
	toprint.append("{0}_{1}".format(sync.chr,sync.pos))
	for a in alcount:
		cov=a[0]+a[1]
		toprint.append(str(cov))
		toprint.append(str(a[1]))
	print "\t".join(toprint)
