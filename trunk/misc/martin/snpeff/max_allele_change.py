import sys,select
import collections
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from optparse import OptionParser,OptionGroup
import copy
import time
import datetime
import math
import random


#Author: Martin Kapun
#version: 2.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python max_allele_change.py -i candidates.sync -p 1,5,9+2,6,10+3,4,8 -s 1,2 -c 10 -m 500 -n 10 -o output
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/ 
3)      this script calculates the allele frequency change of the rising allele summed over all replicates (i.e the sum of all counts for each allele within one timepoint). The results are printed after the cmh p-values
	""") 


parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-p", "--pops", dest="pops", help="define 2 or 3 populations for which you want to calculate the allele frequency change. populations are separated by a \",\"(e.g. 1,5,9+2,6,10+3,4,8)")
parser.add_option("-c", "--count", dest="count", help="minor allele count")
parser.add_option("-m", "--maxcov", dest="maxcov", help="maximum coverage")
parser.add_option("-n", "--mincov", dest="mincov", help="minimum coverage")
parser.add_option("-o", "--out", dest="out", help="output file")
parser.add_option("-s", "--se", dest="se", help="start and end")


parser.add_option_group(group)
(options, args) = parser.parse_args()

def max_allele_change(inp,cp):
	h,h1,h2,h3,h4={},collections.defaultdict(lambda:0),collections.defaultdict(lambda:0),collections.defaultdict(lambda:0),collections.defaultdict(lambda:0)
	cov1,cov2=0,0	
	for r in range(0,len(cp)):
		#print cp[r] 
		cov1+=cov(inp,cp[r][0])
		h1["A"]+=int(inp.split()[int(cp[r][0])+2].split(":")[0])
		h1["T"]+=int(inp.split()[int(cp[r][0])+2].split(":")[1])
		h1["C"]+=int(inp.split()[int(cp[r][0])+2].split(":")[2])
		h1["G"]+=int(inp.split()[int(cp[r][0])+2].split(":")[3])	
		cov2+=cov(inp,cp[r][1])
		h2["A"]+=int(inp.split()[int(cp[r][1])+2].split(":")[0])
		h2["T"]+=int(inp.split()[int(cp[r][1])+2].split(":")[1])
		h2["C"]+=int(inp.split()[int(cp[r][1])+2].split(":")[2])
		h2["G"]+=int(inp.split()[int(cp[r][1])+2].split(":")[3])
	for k,v in h1.items():
		h3[k]=v/float(cov1)
	for k,v in h2.items():
		h4[k]=v/float(cov2)	
	for k,v in h3.items():
		h[k]=h4[k]-v
	if len(set(h.values()))==1.0:
		return "na"
	else:
		maxallele=max(h.values())
		return maxallele
	
def SNP_test(inp,populations,max_coverage,min_coverage,mincount):
	minus,maxcov,mincov=0,0,0
	if l.rstrip()!="":
		for i in range(0,len(populations)): # test if "-" is in sync and coverages for all populations below threshold
			if "-" in inp.split()[int(populations[i])+2] or inp.split()[int(populations[i])+2]=="0:0:0:0:0:0":
				minus=1
				#print "-"
			else:
				if cov(inp,int(populations[i]))>int(max_coverage):
					maxcov=1
					#print "maxcov"
				elif cov(inp,int(populations[i]))<int(min_coverage):
					mincov=1
					#print "mincov"
		if minus!=1 and maxcov!=1 and mincov!=1:
			if minallele_count(inp,populations)!="fixed" and minallele_count(inp,populations)>=int(mincount):
				
				return 1
			else:
				return 0
		else:
			return 0
	return 0

def cov(inp,i): ##returns coverage of population i
	h=collections.defaultdict(lambda:0)
	h["A"]+=int(inp.split()[int(i)+2].split(":")[0])
	h["T"]+=int(inp.split()[int(i)+2].split(":")[1])
	h["C"]+=int(inp.split()[int(i)+2].split(":")[2])
	h["G"]+=int(inp.split()[int(i)+2].split(":")[3])
	cov=sum(h.values())
	if cov!=0:
		return cov
	else:
		return "na"

def minallele_count(inp,pops): ### test wether SNP is "fixed" or below minimum count threshold
	h,h2,h3,mincount,minfreq=collections.defaultdict(lambda:0),[],[],0,0
	for i in range(0,len(pops)):
		h["A"]+=int(inp.split()[int(pops[i])+2].split(":")[0])
		h["T"]+=int(inp.split()[int(pops[i])+2].split(":")[1])
		h["C"]+=int(inp.split()[int(pops[i])+2].split(":")[2])
		h["G"]+=int(inp.split()[int(pops[i])+2].split(":")[3])
	cov=sum(h.values())
	for k,v in h.items():
		if v!=0:
			h2.append(v)
			
	h3=sorted(h2)	
	if len(h3)>2:
		hlim=[]
		while len(hlim)<2:
			hlim.append(h3.pop())
		mincount=min(hlim)
		return mincount
	if len(h3)<2:
		return "fixed"
	else: 
		mincount=min(h3)
		return mincount
	
pops=str(options.pops).split("+")
start=int(str(options.se).split(",")[0])-1
end=int(str(options.se).split(",")[1])-1
o,allpops,comptime=[],[],[]
for line in pops:
	o.append(line.split(","))

genepops=zip(*o)
comprep=[genepops[start],genepops[end]] # replicates at different timepoints
sign=["greater","less"]
signcode={"greater":1,"less":0}
for l in comprep:	
	for i in l:
		allpops.append(i)
for l in pops:
	comptime.append((l.split(",")[start],l.split(",")[end]))
c=0
out=open(options.out,"w")
for l in open(options.inp,"r"):
	c+=1
	if SNP_test(l,allpops,options.maxcov,options.mincov,options.count)!=0:
		if max_allele_change(l,comptime)!="na":
			out.write(l.rstrip()+"\t"+str(max_allele_change(l,comptime))+"\n")
	if c%100000==0:
		print str(c)+" positions processed"
		