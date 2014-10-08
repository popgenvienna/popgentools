import sys
import collections
from optparse import OptionParser, OptionGroup
import copy
import math
import random
from Haplotype import *




#Author: Dr. Robert Kofler
parser = OptionParser()
parser.add_option("--haplotypes", dest="haplotypes", help="the haplotype file")
parser.add_option("--N", dest="n",help="the final population number")
parser.add_option("--output",dest="output", help="the output file")
(options, args) = parser.parse_args()
ofh=open(options.output,"w")

# Determine sample size etc
censussize=int(options.n) # targetsamplesize 
inputsamplesize=HaplotypeIO.haplotypeCount(options.haplotypes) #actualsample size
if(inputsamplesize!=4):
	raise Error("must be two");

for line in open(options.haplotypes):
	line=line.rstrip()
	p=HaplotypeIO.parseLine(line)
	haps=p.haplotypes
	novelhaps=[]
	for i in range(0,censussize):
		novelhaps.append(haps[0])
		novelhaps.append(haps[2])
	# Create novel population from the subsampled loci
	novelp=PopulationHaplotype(p.chrom, p.position, p.refchar, p.major, p.minor, novelhaps) # (self,chrom,position,refchar,major,minor,haplotype):
	topr=HaplotypeIO.formatEntry(novelp)
	ofh.write(topr+"\n")
ofh.close()

