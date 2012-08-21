
import sys
import collections
import modules.RCMH
from modules.CMH import CMHReader
from optparse import OptionParser, OptionGroup
import math


#Author: Martin Kapun


#########################################################   HELP   #########################################################################
#print
usage="python %prog --input full.cmh --all 1,2,3,4,5,6,7,8,9,10 --start 1,2,3 --end 8,9,10 > full.ac"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'''				

H E L P :
_________

Description:
This script returns the  SNP-wise allelel frequency change (AFC) between two time-points (--start and --end) averaged across all replicates (number of populations in --start and --end) and the cumulative coverage (COV) across all populations (in --start and --end). The parameter --all is used to defined the populations used for SNP detection. The two values AFC and COV are appended to the last column of each SNP in the output. Only the two most common alleles (across all populations) will be counted. Per default the script is expecting a sync file as an input. However, it also handles CMH outputs (with a P-value at the end). Then, you need to put --type CMH

''') 
########## 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--all", dest="all", help="all populations used for SNP calling")
parser.add_option("--start", dest="start", help="start population(s) separated by a ',' ")
parser.add_option("--end", dest="end", help="end population(s) separated by a ',' ")
parser.add_option("--type",dest="type",help="file type: cmh output or sync file",default="sync")

parser.add_option_group(group)
(options, args) = parser.parse_args()


def AFC(st,en,allele):
	ac1,ac2,co1,co2,cot=0.0,0.0,0.0,0.0
	for pop in st:
		cot+=pop.cov
		co1+=pop.cov
		ac1+=pop.counth[allele]
	for pop in en:
		cot+=pop.cov
		co2+=pop.cov
		ac2+=pop.counth[allele]
	af=	ac2/co2-ac1/co1
	return af,cot


start=map(int,options.start.split(","))
end=map(int,options.end.split(","))
allpops=map(int,options.all.split(","))


if options.type=="sync":
	filehandle=SyncReader(options.input)
else:
	filehandle=CMHReader(options.input)

for sync in filehandle:
	
	# extract the correct populations
	pops=sync.subpopulations(allpops)
	s=sync.subpopulations(start)
	e=sync.subpopulations(end)
	
	# obtain counts for the two major alleles
	(alcount,ma,mi)=modules.RCMH.Utility.getMajorAlleleCount(pops)
	if AFC(s,e,mi)[0]>0:
		print str(sync)+"\t"+"\t".join(map(str,AFC(s,e,mi)))
	else: 
		print str(sync)+"\t"+"\t".join(map(str,AFC(s,e,ma)))
	
