import sys
import collections
from modules.CMH import PopIO
from optparse import OptionParser, OptionGroup
import copy

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python fixed_SNP.py -i candidates_Base-F15.cmh -p 4,5,6 -t 2
2)	Find fixed SNPs
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-i", "--inp", dest="inp", help="input = cmh file")
parser.add_option("-p", "--pops", dest="pops", help="define replicates and populations: separate populations with a \",\", e.g.8,9,10 ")
parser.add_option("-t", "--th", dest="th", help="minimum allele count ")
parser.add_option("--test",action="store_true", dest="test",help="run the doctest")

parser.add_option_group(group)
(options, args) = parser.parse_args()

def is_fixed(a,pops,th):
	# filter empty lines
	"""
	"""
	for i in pops:
		activepop=a[int(i)-1]
		alcount=activepop.count_alleles(th)
		if alcount>1:
			return 0
	return 1
	

pops=options.pops.split(",")	
for l in open(str(options.inp),"r"):
	if l.rstrip()=="":
		continue
	p=PopIO.parse_cmhline(l)
	isfixed=is_fixed(p.populations,pops,int(options.th))	
	if isfixed>0:
		print l.rstrip()					
	
	