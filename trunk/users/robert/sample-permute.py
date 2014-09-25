import sys
import collections
from optparse import OptionParser, OptionGroup
import copy
import math
import random
from  syncIOpg import *

def get_meshed(populations):
	c=[0,0,0,0]
	for p in populations:
		for i in range(0,4):
			c[i]+=p[i]
	meshed=[]
	for i,a in enumerate(('A','T','C','G')):
		tc=int(c[i])
		[meshed.append(a) for i in range(0,tc)]
	random.shuffle(meshed)
	return meshed
	

def ispolymorphic(populations):
	c=[0,0,0,0]
	for p in populations:
		for i in range(0,4):
			c[i]+=p[i]
	countalleles=0
	for co in c:
		if co>0:
			countalleles+=1
	polymorphic=False
	if countalleles >1:
		polymorphic=True
	return polymorphic

def getcounts(meshed,targetcov):
		counts=[0,0,0,0,0,0] # ATCGNdel
		indexarray={'A':0,'T':1,'C':2,'G':3}
		for i in range(0,targetcov):
			a=meshed.pop()
			index=indexarray[a]
			counts[index]+=1
		return counts
	


	

#Author: Dr. Robert Kofler
parser = OptionParser()
parser.add_option("--input", dest="input", help="the input file")
(options, args) = parser.parse_args()



for chr, pos, t,pops in SyncReaderATCG(options.input):
	if not ispolymorphic(pops):
		continue
	meshedalleles=get_meshed(pops)
	nuevopops=None
	nuevopops=[]
	for p in pops:
		targetcov=int(sum(p[0:4]))
		counts=getcounts(meshedalleles,targetcov)
		nuevopops.append(counts)
	topr=SyncWriterATCG.format(chr,pos,t,nuevopops)
	print topr


	


