import sys
import collections
from optparse import OptionParser,OptionGroup
import copy
from modules.LDUtil import LDIO,Util
import numpy
from modules.CMH import CMHReader

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python LD_dist.py -i candidates.cmh -s full.cmh -b 10 -m 1000 -o out.ld
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script takes the set of SNPs in the proximity of candidates SNPs and spits out the average p-values for SNPs in bins until a maximum distance. E.g. -m 500 and -b 20 will result in bins of 50bp length from -500 to +500bp around the candidate.
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-c", "--candidates", dest="cand", help="*.sync or cmh output file")
parser.add_option("-s", "--snps", dest="snps", help="cmh output with all SNPs")
parser.add_option("-b", "--bins", dest="bins", help="number of bins")
parser.add_option("-m", "--maxdist", dest="maxdist", help="maximum distance to candidate.")
parser.add_option("-o", "--out", dest="out", help="outputfile for boxplot")
parser.add_option("--measure", dest="measure", help="What should be calculated, median (median), geometric mean (gm)")
parser.add_option_group(group)
(options, args) = parser.parse_args()

# 1: Lade die candidaten SNPs
chrh,candl = LDIO.read_candidatehash(options.cand, int(options.maxdist) )

# 2: Iteriere uber alle SNPs
for snp in CMHReader(options.snps):
	chr=snp.chr
	pos=snp.pos
	candidates=chrh[chr][pos]
	for cand in candidates:
		cand.appendSNP(snp)

# 3: itererier ueber alle canidaten
ofh=open(options.out,"w")
for cand in candl:
	bins=cand.distributeToBins(int(options.bins))
	toprint=[]
	toprint.append(cand.chr)
	toprint.append(cand.pos)
	toprint.append(cand.pvalue)
	for t in bins:
		value=0
		pvalar=[snp.pvalue for snp in t]
		if options.measure == "median":
			value=Util.median(pvalar)
		elif options.measure=="gm":
			value=Util.geometricmean(pvalar)
		else:
			raise Exception("Unknown option for --measure: "+options.measure)
	
		if value is None:
			value="na"
		toprint.append(value)
	
	s="\t".join(map(str,toprint))
	ofh.write(s+"\n")

