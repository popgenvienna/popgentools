#!/usr/bin/env python
import sys
import collections
from optparse import OptionParser, OptionGroup
import copy
import math
import random
from Haplotype import *
import gzip

def createRandomIndices(fullSize,subSize):
	toDraw=list(range(0,fullSize))
	toret=[]
	for i in range(0,subSize):
		index=int(random.random()*len(toDraw))
		toret.append(toDraw.pop(index))
	return toret

#Author: Dr. Robert Kofler
parser = OptionParser()
parser.add_option("--haplotypes", dest="haplotypes", help="the haplotype file")
parser.add_option("--2Ne", dest="samplesize",help="the number of chromosomes (haplotypes) to sample")
parser.add_option("--output",dest="output", help="the output file")
(options, args) = parser.parse_args()

f = gzip.open('file.txt.gz', 'wb')
outputfile=options.output
if(not outputfile.endswith(".gz")):
	outputfile+=".gz"
ofh=gzip.open(outputfile, 'wb')
# Determine sample size etc
subsamplesize=int(options.samplesize)
fullsamplesize=HaplotypeIO.haplotypeCount(options.haplotypes)
if(subsamplesize>fullsamplesize):
	raise Exception("subsamplesize nees to be smaller than the number of populations")
	
randindex=createRandomIndices(fullsamplesize,subsamplesize)
print subsamplesize
print fullsamplesize
print randindex

fh=None
file=options.haplotype
if(file.endswith(".gz")):
	fh=gzip.open(file,mode='rb')
else:
	fh=open(file)

for line in fh:
	line=line.rstrip()
	p=HaplotypeIO.parseLine(line)
	haps=p.haplotypes
	novelhaps=[haps[i] for i in randindex]
	
	# Create novel population from the subsampled loci
	novelp=PopulationHaplotype(p.chrom, p.position, p.refchar, p.major, p.minor, novelhaps) # (self,chrom,position,refchar,major,minor,haplotype):
	# Print only polymorphic SNPs
	if(novelp.isPolymorphic):	
		topr=HaplotypeIO.formatEntry(novelp)
		ofh.write(topr+"\n")
ofh.close()
	
