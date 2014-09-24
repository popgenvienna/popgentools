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
	

def getmajmin(freqs,nuevofreqs):
	mii=None
	mai=None
	for i,f in enumerate(freqs):
		if f>0.5:
			mai=i
		elif f >0:
			mii=i
	if mii is None or mai is None:
		return (0,0,0,0)
	toret=(freqs[mii],freqs[mai],nuevofreqs[mii],nuevofreqs[mai])
	return toret


def is_mc_polymorphic(populations,mincount):
	c=[0,0,0,0]
	for p in populations:
		for i in range(0,4):
			c[i]+=p[i]
	cov=float(sum(c))
	countalleles=0
	for co in c:
		if co>=mincount:
			countalleles+=1
	polymorphic=0
	if countalleles >1:
		polymorphic=1

	return polymorphic

#Author: Dr. Robert Kofler
parser = OptionParser()
parser.add_option("--input", dest="input", help="the input file")
parser.add_option("--min-count",dest="mincount",help="the min count")
(options, args) = parser.parse_args()
mincount=int(options.mincount)




for chr, pos, t,pops in SyncReaderATCG(options.input):
	polymorphic,freqs=getfreqs(pops) # A,T,C,G
	if polymorphic:
		ori_mc=is_mc_polymorphic(pops,mincount)
		nuevopops=[]
		for p in pops:
			targetcov=int(sum(p[0:4]))
			counts=getcounts(freqs,targetcov)
			nuevopops.append(counts)
		nuevo_mc=is_mc_polymorphic(nuevopops,mincount)
		tp,nuevofreqs=getfreqs(nuevopops)
		mamif=getmajmin(freqs,nuevofreqs)
		print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(mamif[0],mamif[1],mamif[2],mamif[3],ori_mc,nuevo_mc)



