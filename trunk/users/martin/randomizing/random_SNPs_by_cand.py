import sys
import collections
from optparse import OptionParser,OptionGroup


#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage= "python %prog --cand cand.sync --full full.sync --distance 100000 --out rand_cand.sync"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,"""

H E L P:
____________

This script randomly picks SNPs from an input file (--full; with Chromosome and position in the first two columns) based on the number and poistion of SNPs in the candidate input file (--cand). The script SNPs, which are in an approximate distance (--distance) upstreams of the candidate.

""") 


parser.add_option("-c", "--cand", dest="c", help="any file with Chromosme and positions in the first two columns")
parser.add_option("-f", "--full", dest="f", help="any file with Chromosme and positions in the first two columns")
parser.add_option("-o", "--out", dest="o", help="outfile")
parser.add_option("-d", "--distance", dest="d", help="minimum distance of random SNP to candidate SNP")

parser.add_option_group(group)
(options, args) = parser.parse_args()

#########################################################   CODE   #########################################################################

def positionhash(x):
	poshash=collections.defaultdict(lambda:collections.defaultdict(lambda:""))
	for l in open(x,"r"):
		a=l.split()
		pos=int(a[1])
		chrom=a[0]
		poshash[chrom][pos]=l
	return poshash

cand=positionhash(options.c)
full=positionhash(options.f)
distance=int(options.d)	 
out=open(options.o,"w")
count=0

for chrom,value in cand.items():
	for pos in value.keys():
		distance=int(options.d)	
		count+=1
		if pos-distance<=0:
			while pos+distance not in full[chrom]:
				distance+=1
				#print distance
			out.write(full[chrom][pos+distance])
		else:
			while pos-distance not in full[chrom]:
				distance-=1
				#print distance
			out.write(full[chrom][pos-distance])
