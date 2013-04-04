import sys
import collections
import math
import re
import random
from optparse import OptionParser, OptionGroup



class RHappicker:
	def __init__(self,dgrpsize):
		self.__templist=range(0,dgrpsize)
		self.__dgrpsize=dgrpsize
		
	def __newlist(self):
		return list(self.__templist)
	
	def get_random(self,popsize,homozygous):
		target=popsize*2
		if(homozygous and popsize<self.__dgrpsize):
			raise ValueError("not valid population size")
		if((not homozygous) and popsize<(self.__dgrpsize/2)):
			raise ValueErro("not valid population size")
		toret=[]
		todraw=self.__newlist()
		
		while(len(toret)<target):
			if(len(todraw)==0):
				todraw=self.__newlist()
			r=int(random.random()*len(todraw))
			e=todraw.pop(r)
			if(homozygous):
				toret.append(e)
			toret.append(e)
		return toret


def pick_genotypes(genotypelist,randhaplotypes):
	toret=[]
	for r in  randhaplotypes:
		toret.append(genotypelist[r])
	return toret

def format_haplotypes(genotypelist):
	gsize=len(genotypelist)
	toformat=[]
	for i in range(0,gsize,2):
		toformat.append(genotypelist[i]+genotypelist[i+1])
	toret=" ".join(toformat)
	return toret

def get_allele_hash(genotypelist):
	toret=collections.defaultdict(lambda:0)
	for g in genotypelist:
		toret[g]+=1
	return toret

def format_majorminor(allelehash):
	ts=allelehash.items()
	ts=sorted(ts,key=lambda a: a[1])
	toret=ts[1][0]+"/"+ts[0][0]
	return toret
	
	
### INPUT
# X	207	C	NCNNCNNCNNCNNCCCCCNNNNNNNNNNCNCCNNCNNNNNNCNNNCNNNCNNNNCNNCCCNNNNCANNNCNNNNNCNNNNNCNCNNNCNNNNNNNNNNNNNNNNNCCCNNNNNNNCNNCNNNCNNNNCCNNNNCNCNCNNNCNNNNCNCNCNNCCNNC
# X	208	A	NANNANNANNANNAAGANNNNNNNNNNNANAANNANNNNNNANNNANANAANANANNAAANNNNAANNNANNNNNANANNNANANNNANNNNNNNNANNNNNNNNAAANNNNNNNANNANNNANNNNAANNNNNNANANNNANNNNANANANNNANNA
# X	221	G	NGNNANNGNNGGNGGGGNNGNNGNNNNNGNGGNNGGNNNTNGNNNGNGNNGGGNGNNNGGNNNNGGNNNGNNNNNGGGNNNGNGNGNGNNNNNGNNGGNNNNNNNGGGNNNNNNGGNNGNNNNNNNNGGNNNNGNGNGNNNGGNNNGGGNGNNGGNNG


### OUTPUT
# 2L      686891    G      G/A    AG GG GG GG GG
# 2L      936681    A      A/G    GG AA AA AA AA
# 2L      861026    T      A/T    TT AT AA AA TT
# 2L      966618    C      T/C    TC TC TT TT TT
# 2R      134298    A      A/C    AC AC CC CC CC


# Author: Dr. Robert Kofler
#
# Description
# Creates an 'N' -allele frequency spectrum for dgrp haplotype table
# This is useful to decide which SNPs may be used and which ones not
# Note: the script treats everything not ATCG as N (including M and other nonsenese)!

dgrpsize=158
parser = OptionParser()
parser.add_option("--input",dest="input",help="A file containing a haplotype table")
parser.add_option("--population-size",dest="popsize", help="The size of the population to build")
parser.add_option("--homozygot",dest="homozygote",action="store_true",help="Should a homozygous population be created; default=heterozygous")
(options, args) = parser.parse_args()
popsize=int(options.popsize)
homozygous=bool(options.homozygote)

rpicker=RHappicker(dgrpsize)
randgenotypes=rpicker.get_random(popsize,homozygous)


for line in open(options.input):
	line=line.rstrip('\n')
	a=line.split('\t')
	genotypes=pick_genotypes(list(a[3]),randgenotypes)
	allelehash=get_allele_hash(genotypes)
	if(len(allelehash)<2):
		continue
	
	topr=[]
	topr.append(a[0])
	topr.append(a[1])
	topr.append(a[2])
	topr.append(format_majorminor(allelehash))
	topr.append(format_haplotypes(genotypes))
	tostr="\t".join(topr)
	print tostr
	
	