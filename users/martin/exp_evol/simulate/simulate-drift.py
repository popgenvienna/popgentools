#!/usr/bin/env python
import sys
import collections
import modules.RCMH
from modules.CMH import SyncReader
from optparse import OptionParser, OptionGroup


#Author: Robert Kofler


#########################################################   HELP   #########################################################################
#print
usage="python %prog --input full.sync --output full_random.sync --min-count 2 --min-coverage 20 --maxcoverage 500 --base-index 1,2,3 --derived-index 4,5,6,7,8,9 --derived-timestamps 15,15,15,37,37,37 --populations-size 200"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'''				

H E L P :
_________

Description:

This script performs Wright Fisher sampling to simulate allele frequencies changes only due to genetic drift. This can be done in parallel for multiple populations, which is useful if done for multiple replicates. Therefore, the columns of the replicates in the sync file must be defined with --base-index. Then the allele frequencies of these populations will be used as the starting frequencies for the simulations. The parameter --derived-timestamps defines the number of generations of simulations for each replicate. In case you have three replicates and want to simulate 15 generations you need to put "--derived-timestamps 15,15,15". The parameter --derived-index defines the columns of the populations in the input file which are corresponding to the simulated generations. The simulations will be subsampled to the same coverage as in these populations. --population-size defines the effective populations size used for the simulations. 
The output is a synchronized file containig the Basepopulations defined with --base-index and the simulated population defined with --derived-timestamps.

	''') 
########## 
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
parser.add_option_group(group)
(options, args) = parser.parse_args()

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
	
	# create new single populations and a total new population line
	singlepops=modules.RCMH.Utility.getSinglePops(alcount+simulated,majora,minora)
	newp=modules.CMH.PopLine(sl.chr,sl.pos,sl.refc,singlepops)


	ofw.write(str(newp)+"\n")


	
