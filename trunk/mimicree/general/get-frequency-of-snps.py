import sys
import collections
import math
from optparse import OptionParser, OptionGroup
import gzip
from mimicreeHapIO import MimicreeHapReader

def loadsnps(file):
	toret={}
	for l in open(file):
		l=l.rstrip()
		a=l.split("\t")
		key=a[0]+":"+a[1]
		nonselected=a[2]
		toret[key]=nonselected
	return toret


parser = OptionParser()
parser.add_option("--haplotype",dest="haplotype",help="the haplotype for which the frequency should be displayed")
parser.add_option("--snps",dest="snps",help="The sync or cmh file")
(options, args) = parser.parse_args()

snps=loadsnps(options.snps)

for mh in MimicreeHapReader(options.haplotype):
	key=mh.chr+":"+str(mh.pos)
	if key in snps:
		nsa=snps[key]
		freq_not_selected=mh.freqAllele(nsa)
		freq_selected=1.0-freq_not_selected
		print "{0}\t{1}".format(key,freq_selected)
	