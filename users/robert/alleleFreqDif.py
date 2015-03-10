import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from  syncIOpg import SyncReaderMajMin
import math




def parseComparision(compstr):
	topa=[compstr]
	if "," in compstr:
		topa=compstr.split(",")
	toret=[]
	for t in topa:
		a,b=t.split("-")
		toret.append((int(a),int(b)))
	return toret

def transformFreq(totrans):
	cov=(totrans[0]+totrans[1])
	if cov==0:
		return (0.0,0.0)
	cov=float(cov)	
	af=float(totrans[0])/cov
	bf=float(totrans[1])/cov
	return (af,bf)
	

parser = OptionParser()

parser.add_option("--input",dest="input",help="the input file as sync")
parser.add_option("--min-count",dest="mincount",help="the minimum frequency")
parser.add_option("--min-cov",dest="mincov",help="the minimum coverage")
parser.add_option("--to-compare",dest="tocomp",help="samples to compar 1-2,3-4,5-6,7-8")

(options, args) = parser.parse_args()
toCompare=parseComparision(options.tocomp)
mincount=int(options.mincount)
mincov=int(options.mincov)

fdH=collections.defaultdict(lambda:0)
nH=collections.defaultdict(lambda:0)

for chr,pos,temp,mami in SyncReaderMajMin(options.input):
	for toco in toCompare:
		aindex=toco[0]
		bindex=toco[1]
		a=mami[aindex-1]
		b=mami[bindex-1]
		af=transformFreq(a)
		bf=transformFreq(b)
		if(a[0]>=mincount and a[1]>=mincount and b[0]>=mincount and b[1]>=mincount and a[0]+a[1]>=mincov and b[0]+b[1]>=mincov):
			amf=af[0] # major frequency of first
			bmf=bf[0] # major frequency of second
			fd=math.fabs(amf-bmf)
			fdH[toco]+=fd
			nH[toco]+=1 # counts for comparision



print "to compare {0}".format(toCompare)
print "min coverage {0}".format(mincov)
print "min count {0}".format(mincount)
for k in nH.keys():
	n=nH[k]
	fdsum=fdH[k]
	avdiv=fdsum/n

	print "{0}\t{1}\t{2}\t{3}".format(k,n,fdsum,avdiv)
	
		






