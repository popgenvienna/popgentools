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

sH=collections.defaultdict(lambda:{"a":0,"b":0})
ssH=collections.defaultdict(lambda:{"a":0,"b":0})
ssxyH=collections.defaultdict(lambda:0)


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
			nH[toco]+=1 # counts for comparision
			sH[toco]["a"]+=amf
			sH[toco]["b"]+=bmf
			ssH[toco]["a"]+=amf*amf
			ssH[toco]["b"]+=bmf*bmf
			ssxyH[toco]+=amf*bmf

for k in nH.keys():
	aindex=k[0]
	bindex=k[1]
	n=nH[k]
	sumx=sH[k]["a"]
	sumy=sH[k]["b"]
	sumxx=ssH[k]["a"]
	sumyy=ssH[k]["b"]
	sumxy=ssxyH[k]
	ssxy=sumxy-((sumx*sumy)/n)
	ssx=sumxx-((sumx*sumx)/n)
	ssy=sumyy-((sumy*sumy)/n)
	r=ssxy/math.sqrt(ssx*ssy)
	print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(k,n,ssx,ssy,ssxy,r)
	
		






