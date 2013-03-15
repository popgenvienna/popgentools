import sys
import collections
from optparse import OptionParser, OptionGroup
import copy
import math
import random
from Haplotype import *


def createRandomIndices(inputSize,hapCopies):
	toDraw=range(0,inputSize)
	toret=[]
	while(len(toDraw)>0):
		idx=int(random.random()*len(toDraw))
		drawn=toDraw.pop(idx)
		for i in range(0,hapCopies):
			toret.append(drawn)
	return toret


#Author: Dr. Robert Kofler
parser = OptionParser()
parser.add_option("--haplotypes", dest="haplotypes", help="the haplotype file")
parser.add_option("--2Ne", dest="samplesize",help="the final population number; The number of chromosomes in the input has to be less than the number of chromosomes in the output")
parser.add_option("--output",dest="output", help="the output file")
(options, args) = parser.parse_args()
ofh=open(options.output,"w")

# Determine sample size etc
boostsamplesize=int(options.samplesize) # targetsamplesize 
inputsamplesize=HaplotypeIO.haplotypeCount(options.haplotypes) #actualsample size
if(boostsamplesize<=inputsamplesize):
	raise Exception("boostsamplesize needs to be larger than the number of haplotypes in the input" )
hapCopies=int(float(boostsamplesize)/float(inputsamplesize))
if(hapCopies%2!=0 or hapCopies*inputsamplesize != boostsamplesize):
	raise Exception("boostsamplesize needs to be a twofold of the inputsample size")
	

randindex=createRandomIndices(inputsamplesize,hapCopies)
print boostsamplesize
print inputsamplesize
print randindex

for line in open(options.haplotypes):
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

