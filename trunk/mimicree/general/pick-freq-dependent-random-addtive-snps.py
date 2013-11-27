import sys
import collections
import math
import re
import random
from optparse import OptionParser, OptionGroup


# Author: Dr. Robert Kofler
# Description: Picks selected SNPs from Andreas simulated natural population, fastcoalsim
# Update: the major/minor distinction between alleles has now been deprecated; we are now working with ancestral(=mostly major) and derived(=mostly minor)

class SelectedCandidate:
	def __init__(self,chr,pos,w11,freq):
		self.chr=chr
		self.pos=pos
		self.w11=w11
		self.freq=freq

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

	derivedfreq=float(derivedcount)/float(totcount)	
	return (chr,pos,ancestral,derived,derivedfreq)


parser = OptionParser()
parser.add_option("--input",dest="input",help="A file containing the dgrp haplotypes (MimicrEE input)")
parser.add_option("-n","--number-selected",dest="numsel",help="Number of selected SNPs")
parser.add_option("-e","--heterozygous",dest="het",help="Heterozygous effect")
parser.add_option("-s","--selection",dest="selcoef",help="Selection coefficient")
parser.add_option("--min-frequency",dest="minfrequency", help="The minimum frequency")
parser.add_option("-m","--max-frequency",dest="maxfrequency",help="Max. frequency of selected allele")
parser.add_option("--loci-count",dest="locicount",help="Approximate loci count; Necessary for random picking approximation")
(options, args) = parser.parse_args()

apcount=float(options.locicount)/1000.0
numsel=int(options.numsel)
pickfreq=float(numsel)/apcount
maxfreq=float(options.maxfrequency)
minfreq=float(options.minfrequency)

h=float(options.het)
s=float(options.selcoef)

cand=[]
for line in open(options.input):
	if(random.random()>pickfreq):
		continue
		# process only loci that will end up in the  initial set
	line=line.rstrip()
	(chr,pos,ancestral,derived,derivedfreq)=parse_line(line)
	
	
	if(random.random()<0.5):
		# the derived is selected
		if(derivedfreq < maxfreq and derivedfreq >= minfreq):
			cand.append(SelectedCandidate(chr,pos,ancestral,derivedfreq)) # w11=ancestral # w22=selected = derived
	else:
		# the ancestral is selected
		ancfreq=1.0-derivedfreq
		if(ancfreq < maxfreq and ancfreq >= minfreq):
			cand.append(SelectedCandidate(chr,pos,derived,ancfreq)) # w11=derived # w22=ancestral


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
	o.append(str(t.freq))
	tmp="\t".join(o)
	print tmp