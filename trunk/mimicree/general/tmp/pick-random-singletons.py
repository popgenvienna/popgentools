import sys
import collections
import math
import re
import random
from optparse import OptionParser, OptionGroup
import gzip


# Author: Dr. Robert Kofler
# Description: Picks selected SNPs from Andreas simulated natural population, fastcoalsim
# Update: the major/minor distinction between alleles has now been deprecated; we are now working with ancestral(=mostly major) and derived(=mostly minor)

class SelectedCandidate:
	def __init__(self,chr,pos,w11):
		self.chr=chr
		self.pos=pos
		self.w11=w11

def parse_line(line):
	"""
	2L      686891    G      G/A    AG GG GG GG GG
	2L      936681    A      A/G    GG AA AA AA AA
	2L      861026    T      A/T    TT AT AA AA TT
	SCHEISSE:
	2L 	230000	  T	 A/T	TT TT TT TT TT
	"""
	a=line.split("\t")
	chr=a[0]
	pos=a[1]
	ancestral,derived=a[3].split("/")
	b=a[4].split(" ")
	nucs=[]
	for i in b:
		nucs.extend(list(i))
	derivedcount=0
	totcount=0
	for n in nucs:
		totcount+=1
		if(n==derived):
			derivedcount+=1
	ancestralcount=totcount-derivedcount
	
	# define major and minor allele
	maja,mina=(ancestral,derived)
	majc,minc=(ancestralcount,derivedcount)
	if(derivedcount>ancestralcount):
		maja,mina=(derived,ancestral)
		majc,minc=(derivedcount,ancestralcount)

	# (chr,pos,minorallele,majorallele,minorcount)=parse_line(line)
	return (chr,pos,mina,maja,minc)


parser = OptionParser()
parser.add_option("--input",dest="input",help="A file containing the dgrp haplotypes (MimicrEE input)")
parser.add_option("-n","--number-selected",dest="numsel",help="Number of selected SNPs")
parser.add_option("-e","--heterozygous",dest="het",help="Heterozygous effect")
parser.add_option("-s","--selection",dest="selcoef",help="Selection coefficient")
parser.add_option("--singleton-definition",dest="singledef",help="Exact count of suitable SNPs")
parser.add_option("--loci-count",dest="locicount",help="Approximate loci count; Necessary for random picking approximation")
(options, args) = parser.parse_args()

apcount=float(options.locicount)/20.0
numsel=int(options.numsel)
pickfreq=float(numsel)/apcount
singledef=float(options.singledef)

h=float(options.het)
s=float(options.selcoef)

cand=[]
file=options.input
fh=None
if(file.endswith(".gz")):
	fh=gzip.open(file,mode='rb')
else:
	fh=open(file)

for line in fh:
	if(random.random()>pickfreq):
		continue
		# process only loci that will end up in the  initial set
	line=line.rstrip()
	(chr,pos,minorallele,majorallele,minorcount)=parse_line(line)
	
	if minorcount!=singledef:
		continue

	cand.append(SelectedCandidate(chr,pos,majorallele))
	# provide the not-selected: ie the major
	# w11=major w22=minor=selected




if(len(cand) < numsel):
	raise ValueError("not enough SNPs initially picked; please provide a more realistic estimate of the amount of loci present (~factor 10) "+len(cand))

toprint=[]
while(len(toprint)<numsel):
	r=int(random.random()*len(cand))
	toprint.append(cand.pop(r))
	

for t in toprint:
	o=[]
	o.append(t.chr)
	o.append(t.pos)
	o.append(t.w11)
	o.append(str(s))
	o.append(str(h))
	tmp="\t".join(o)
	print tmp