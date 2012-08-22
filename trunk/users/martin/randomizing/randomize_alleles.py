import sys
import collections
import modules.RCMH
from modules.CMH import SyncReader
from optparse import OptionParser, OptionGroup
import random


#Author: Martin Kapun


#########################################################   HELP   #########################################################################
usage= "python %prog --input.sync  > randomized_input.sync"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

This script randomly assigns alleles to each population after pooling all populations. This is done by random sampling with replacement. The coverages of each population stays the same.
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A sync file")
parser.add_option_group(group)
(options, args) = parser.parse_args()

synclist=["A","T","C","G","N","Del"]

apops=range(1,len(open(options.input).readline().split()[3:])+1)
#print apops

for sync in SyncReader(options.input):
	
	# extract the correct populations
	pops=sync.subpopulations(apops)
	
	if not modules.RCMH.Utility.coverageValid(pops, 1, 1000000):
		continue
		
	# obtain counts for the two major alleles
	(alcount,majora,minora)=modules.RCMH.Utility.getMajorAlleleCount(pops)
	minc,popcov=0,[]
	for item in pops:
		popcov.append(item.cov)
		minc+=item.counth[minora]
	if minc==0:
		print sync.__str__()
	else: 
		full_a=minc
		full_A=sum(popcov)-full_a
		dictcov=collections.defaultdict(lambda:0)
		syn=[]
		for i in range(sum(popcov)):
			if i < sum(popcov)-full_a:
				dictcov[i]=0
			else:
				dictcov[i]=1
		for pop_coverage in popcov:
			randcov=collections.defaultdict(lambda:0)
			randkeys=random.sample(dictcov.keys(),pop_coverage)
			for k in randkeys:
				randcov[k]=dictcov[k]
				del dictcov[k]
			new_a=sum(randcov.values())
			syncl=[]
			for nuc in synclist:
				if nuc==majora:
					syncl.append(str(int(pop_coverage-new_a)))
				elif nuc==minora:
					syncl.append(str(int(new_a)))
				else:
					syncl.append("0")
			syn.append(":".join(syncl))
		print sync.chr+"\t"+str(sync.pos)+"\t"+sync.refc+"\t"+"\t".join(syn)
			
		
	

	