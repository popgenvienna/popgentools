import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from  syncIOpg import SyncReaderMajMin
import math
import heapq



class FreqContainer:
	def __init__(self,maxsize):
		self.__afd=heapq.heapify([])
		self.__maxsize=maxsize
		self.__minafd=0.0
	
	def add(self,afd):
		if afd< self.__minafd:
			return 0
		
		if len(self.__afd) < self.__maxlen:
			self
		

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
parser.add_option("--min-afd",dest="minafd",help="the minimum allele frequency difference")

(options, args) = parser.parse_args()
toCompare=parseComparision(options.tocomp)
mincount=int(options.mincount)
mincov=int(options.mincov)
minafd=float(options.minafd)


afd=collections.defaultdict(lambda:[])
# 1-3,2-4,1-5,2-6,1-7,2-8,3-5,4-6,3-7,4-8,5-7,6-8
# 1-2,3-4,5-6,7-8,1-3,2-4,1-5,2-6,1-7,2-8,3-5,4-6,3-7,4-8,5-7,6-8
afc=collections.defaultdict(lambda:0)

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
			if fd>minafd:
				afd[toco].append(fd)
			afc[toco]+=1



print "to compare {0}".format(toCompare)
print "min coverage {0}".format(mincov)
print "min count {0}".format(mincount)
for k,v in afc.items() :
	print "@totcounts\t{0}\t{1}".format(k,v)
	
for k in afd.keys():
	h=afd[k]
	for v in h:		
		print "{0}\t{1}".format(k,v)
	
		

#afd_rep = [0.0 for i in range(0,size)]
#afd_xp  = [0.0 for i in range(0,size)]
#reps=set([(1,2),(3,4),(5,6),(7,8)])
#xps=set([(1,3),(2,4),(1,5),(2,6),(1,7),(2,8),(3,5),(4,6),(3,7),(4,8),(5,7),(6,8)])




