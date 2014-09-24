import sys
import collections
from optparse import OptionParser, OptionGroup
import copy
import math
import random
from  syncIO import * 

def getfreqs(populations):
	c=[0,0,0,0]
	for p in populations:
		for i in range(0,4):
			c[i]+=p[i]
	cov=float(sum(c))
	countalleles=0
	for co in c:
		if co>0:
			countalleles+=1
	freqs=[float(i)/cov for i in c]	
	polymorphic=False
	if countalleles >1:
		polymorphic=True

	return polymorphic,freqs

def getcounts(freqs,targetcov):
		tuple=(freqs[0], freqs[0]+freqs[1], freqs[0]+freqs[1]+freqs[2], 1.0)
		counts=[0,0,0,0,0,0] # ATCGNdel
		for c in range(0,targetcov):
			r=random.random()
			index=None
			for i,t in enumerate(tuple):
				if r < t:
					index=i
					break
			counts[index]+=1
		return counts
	


	

#Author: Dr. Robert Kofler
parser = OptionParser()
parser.add_option("--input", dest="input", help="the input file")
(options, args) = parser.parse_args()



for chr, pos, t,pops in SyncReaderATCG(options.input):
	polymorphic,freqs=getfreqs(pops) # A,T,C,G
	nuevopops=None
	if polymorphic:
		nuevopops=[]
		for p in pops:
			targetcov=int(sum(p[0:4]))
			counts=getcounts(freqs,targetcov)
			nuevopops.append(counts)

	else:
		nuevopops=pops
			
	topr=SyncWriterATCG.format(chr,pos,t,nuevopops)
	print topr
	


