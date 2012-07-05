import sys
import collections
import math
import numpy
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version: 2.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: filter_mea.py -i candidates.cmh -f 0.1 -p 1,5,9+2,6,10+3,4,8 -s 1,2 -o out -c 0.1
2)      filter dataset by average minimum allele frequency change (-c) and minimum initial allele frequency (-f)
	""") 
#########################################################   CODE   #########################################################################


parser.add_option("-i", "--inp1", dest="inp1", help="input1: cmh or sync file")
parser.add_option("-f", "--freq", dest="freq", help="minimum initial allele frequeny of allele rising in frequency")
parser.add_option("-p", "--pops", dest="pops", help="define, for which population in the sync you want to produce the minor allele frequency distribution")
parser.add_option("-s", "--se", dest="se", help="define start and end timepoints in analysis, e.g. Base-F15 = 1,2")
parser.add_option("-o", "--out", dest="out", help="out file")
parser.add_option("-c", "--change", dest="change", help="allele-frequency change cutoff")

parser.add_option_group(group)
(options, args) = parser.parse_args()

############################################## functions ############################################################
def freqhash(inp,i):
	h,h2=collections.defaultdict(lambda:0),{}
	h["A"]+=int(inp.split()[i+2].split(":")[0])
	h["T"]+=int(inp.split()[i+2].split(":")[1])
	h["C"]+=int(inp.split()[i+2].split(":")[2])
	h["G"]+=int(inp.split()[i+2].split(":")[3])
	cov=sum(h.values())
	if cov!=0:
		for k,v in h.items():
			h2[k]=float(v)/cov
		#print h2
		return h2
	else:
		return "na"
	
def meanstdv(x): ### calculate mean, median, stdev. standard error : x=datalist
	from math import sqrt
	n,mean,std,se,median = len(x),0,0,0,0
	for a in x:
		mean = mean + a
	mean = mean / float(n)
	for a in x:
		std = std + (a - mean)**2
	std = sqrt(std / float(n-1))
	se= std/sqrt(n)
	median=sorted(x)[len(x)/2]
	return mean,std,se,median

def max_allele_change(inp,cp):
	h,h1,h2={},collections.defaultdict(lambda:0),collections.defaultdict(lambda:0)
	for i in cp[0]:
		#print i
		h1["A"]+=int(inp.split()[int(i)+2].split(":")[0])
		h1["T"]+=int(inp.split()[int(i)+2].split(":")[1])
		h1["C"]+=int(inp.split()[int(i)+2].split(":")[2])
		h1["G"]+=int(inp.split()[int(i)+2].split(":")[3])
	for i in cp[1]:
		#print i
		h2["A"]+=int(inp.split()[int(i)+2].split(":")[0])
		h2["T"]+=int(inp.split()[int(i)+2].split(":")[1])
		h2["C"]+=int(inp.split()[int(i)+2].split(":")[2])
		h2["G"]+=int(inp.split()[int(i)+2].split(":")[3])
	for k,v in h1.items():
		h[k]=h2[k]-v
	if len(set(h.values()))==1.0:
		return "na"
		print "na"
	else:
		maxallele=max(h,key=lambda a: h.get(a))
		return maxallele
		
def allelefreq_diff(inp,change,cp,out,cr,freqth):
	outth=open(out+"_freq_"+str(freqth)+"-change_"+str(change),"w")
	for l in open(inp,"r"):
		diffl,freqs=[],[]
		mea=max_allele_change(l,cr)
		for i in range(0,len(cp)):
			pop1freq=freqhash(l,int(cp[i][0]))[mea]
			pop2freq=freqhash(l,int(cp[i][1]))[mea]
			diffl.append(pop2freq-pop1freq)
			freqs.append(pop1freq)
		if meanstdv(diffl)[0]>=change and min(freqs)>=freqth:
			outth.write(l)
			

	
############################################## code ############################################################	
freq=float(options.freq)
change=float(options.change)
pops=str(options.pops).split("+")
start=int(str(options.se).split(",")[0])-1
end=int(str(options.se).split(",")[1])-1
o,allpops,comptime=[],[],[]
for line in pops:
	o.append(line.split(","))

genepops=zip(*o)
comprep=[genepops[start],genepops[end]] 																						# replicates at different timepoints
for l in comprep:	
	for i in l:
		allpops.append(i)
for l in pops:
	comptime.append((l.split(",")[start],l.split(",")[end])) 
	
allelefreq_diff(options.inp1,change,comptime,options.out,comprep,freq)