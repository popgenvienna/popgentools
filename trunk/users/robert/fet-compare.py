import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from  syncIOpg import SyncReaderMajMin
import math
from fisher import pvalue




def parseComparision(compstr):
	topa=[compstr]
	if "," in compstr:
		topa=compstr.split(",")
	toret=[]
	for t in topa:
		a,b=t.split("-")
		toret.append((int(a),int(b)))
	return toret

	

parser = OptionParser()

parser.add_option("--input",dest="input",help="the input file as sync")
parser.add_option("--min-count",dest="mincount",help="the minimum frequency")
parser.add_option("--min-cov",dest="mincov",help="the minimum coverage")
parser.add_option("--to-compare",dest="tocomp",help="samples to compar 1-2,3-4,5-6,7-8")

(options, args) = parser.parse_args()
toCompare=parseComparision(options.tocomp)
mincount=int(options.mincount)
mincov=int(options.mincov)




for chr,pos,temp,mami in SyncReaderMajMin(options.input):
	for toco in toCompare:
		aindex=toco[0]
		bindex=toco[1]
		a=mami[aindex-1]
		b=mami[bindex-1]
		if(a[0]>=mincount and a[1]>=mincount and b[0]>=mincount and b[1]>=mincount and a[0]+a[1]>=mincov and b[0]+b[1]>=mincov):
			p=pvalue(a[0],a[1],b[0],b[1])
			print toco, p.two_tail



	
		






