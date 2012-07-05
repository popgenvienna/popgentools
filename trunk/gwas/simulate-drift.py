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
Description

Simulate drift using the allele frequencies of a given population
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--output", dest="output", help="A output file")
parser.add_option("--min-count",dest="mincount",help="The minimum allele count")
parser.add_option("--min-coverage",dest="mincoverage",help="The minimum coverage")
parser.add_option("--max-coverage",dest="maxcoverage",help="The maximum coverage")
parser.add_option("--base-index",dest="basecolumn",help="The indices of the base population in the synchronized file,e.g.: 1,2,3")
parser.add_option("--derived-index",dest="derivedcov",help="The indices of the derived populations in the synchronized file. This will only be used for the coverage. E.g.: 4,5,6,7,8,9 or 4,4,4,4,4,4")
parser.add_option("--derived-timestamps",dest="derivedsnap",help="The time stamps at which the allele frequencies of the derived populations should be samples; Has to be in the same order as the derived indices; e.g.: 15,15,23,37,37,37")
parser.add_option("--population-size",dest="popsize",help="Population size in chromosomes (-> provide number of diploids*2)")
parser.add_option("--test",action="store_true", dest="test",help="run the doctest")
parser.add_option_group(group)
(options, args) = parser.parse_args()

# --input /Users/robertkofler/dev/testfiles/base.sync --output /Users/robertkofler/dev/testfiles/simulate.txt --min-count 3 --min-coverage 10 --max-coverage 400 --base-column 1,2,3 --derived-coverages 4,5,6,8,9,10 --derived-snapshots 15,15,23,37,37,37 --population-size 200

mincount	= int(options.mincount)
mincoverage	= int(options.mincoverage)
maxcoverage	= int(options.maxcoverage)
basepopulations	= map(lambda x: int(x),options.basecolumn.split(","))
derivedcolumn	= map(lambda x: int(x),options.derivedcov.split(","))
derivedtimestamp= map(int,options.derivedsnap.split(","))
popsize		= int(options.popsize)
if len(derivedcolumn)!=len(derivedtimestamp):
	raise IOError("Invalid arguments; derived coverages have to have the same length as the derived timestamps")
if len(derivedcolumn)%len(basepopulations)!=0:
	raise IOError("Invalid arguments; derived coverages has to a multiple of the basepopulations")


ofw=open(options.output,"w")
for sl in SyncReader(options.input):
	
	# extract the base populations
	base=sl.subpopulations(basepopulations)
	# test if coverage is valid
	if not modules.RCMH.Utility.coverageValid(base, mincoverage, maxcoverage):
		continue
	# exclude deletions
	delcount=modules.RCMH.Utility.getDeletionCount(base)
	if delcount>=1:
		continue
	# obtain counts for the two major alleles
	(alcount,majora,minora)=modules.RCMH.Utility.getMajorAlleleCount(base)
	# check if SNP
	if not modules.RCMH.Utility.issnp(alcount,mincount):
		continue
	
	# get derived coverages
	derived=sl.subpopulations(derivedcolumn)
	derivedcoverages=map(lambda d:d.cov,derived)
	
	simulated=modules.RCMH.SimulateDrift.multiSampleDrift(popsize,alcount,derivedtimestamp,derivedcoverages)
	#1000,((20,20),(0,20),(20,0)),(20,20,20,40,40,40),(10,20,30,40,50,60) -> [(8, 2), (0, 20), (30, 0), (30, 10), (0, 50), (60, 0)]
	
	# create new single populations and a total new population line
	singlepops=modules.RCMH.Utility.getSinglePops(alcount+simulated,majora,minora)
	newp=modules.CMH.PopLine(sl.chr,sl.pos,sl.refc,singlepops)

	
	# print it!
	ofw.write(str(newp)+"\n")


	
