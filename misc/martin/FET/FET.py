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

1)	usage: python FET_stouffer.py -i candidates.sync -b 100 -p 1,5,9+2,6,10+3,4,8 -s 1,2 -b 100 -t 0.1 -c 10 -m 500 -n 10 -o output
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/ 
	""") 


parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-p", "--pops", dest="pops", help="define populations for which you want to calculate the allele frequency change. replicate populations are separated by a \",\" and different treatments or timeposints are separated by a \"+\" (e.g. 1,5,9+2,6,10+3,4,8)")
parser.add_option("-t", "--th", dest="th", help="minor allele frequency threshold")
parser.add_option("-c", "--count", dest="count", help="minor allele count threshold summed across all populations")
parser.add_option("-m", "--maxcov", dest="maxcov", help="maximum coverage threshold for each population separatehd")
parser.add_option("-n", "--mincov", dest="mincov", help="minimum coverage")
parser.add_option("-o", "--out", dest="out", help="output file")


parser.add_option_group(group)
(options, args) = parser.parse_args()

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
		h["A"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[0])
		h["T"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[1])
		h["C"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[2])
		h["G"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[3])
	cov=sum(h.values())
	for k,v in h.items():
		if v!=0.0:
			h2.append(v)
	h3=sorted(h2)	
	if len(h3)>1:
		hlim=[]
		while len(hlim)<2:
			hlim.append(h3.pop())
		mincount=min(hlim)
		return mincount
	if len(h3)<=1:
		return "fixed"

def counthash(inp,i): ### returns allelefrequencies fo a given population i
	h=collections.defaultdict(lambda:0)
	h["A"]+=int(inp.split()[int(i)+2].split(":")[0])
	h["T"]+=int(inp.split()[int(i)+2].split(":")[1])
	h["C"]+=int(inp.split()[int(i)+2].split(":")[2])
	h["G"]+=int(inp.split()[int(i)+2].split(":")[3])
	return h
	
def SNP_test(inp,populations,max_coverage,min_coverage,mincount,minfreq):
	minus,maxcov,mincov=0,0,0
	if l.rstrip()!="":
		for i in range(0,len(populations)): # test if "-" is in sync and coverages for all populations below threshold
			if "-" in inp.split()[int(populations[i])+2] or "0:0:0:0:0:0" in inp.split()[int(populations[i])+2]:
				minus=1
			else:
				if cov(inp,int(populations[i]))>int(max_coverage):
					print cov(inp,int(populations[i])),int(max_coverage),int(min_coverage)
					maxcov=1
				elif cov(inp,int(populations[i]))<int(min_coverage):
					mincov=1		
		if minus!=1 and maxcov!=1 and mincov!=1:
			if minallele_count(inp,populations)!="fixed" and minallele_count(inp,populations)>=int(mincount):
				return 1
			else:
				return 0
		else:
			return 0
	return 0

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

def hist(x,bins,out,title): ##produce a histogram output: x=data, bins=Bins, out=outputfile, title=Title
	import rpy2.robjects as robjects
	from rpy2.robjects.packages import importr
	r=robjects.r
	graphics = importr('graphics')
	grdevices = importr('grDevices')
	#rang=robjects.vectors.FloatVector([-4.75,-4.74])
	values=robjects.vectors.FloatVector(x)
	grdevices.png(out+'_hist.png')
	graphics.hist(values,breaks=bins,main=title+"\n\n"+str(bins)+" bins",xlab="")
	grdevices.dev_off()


def FET_perm(x):
	r=robjects.r
	popcounts = robjects.IntVector(x)
	popmatrix = robjects.r['matrix'](popcounts, nrow = 2,byrow=False)
	FETres = [0.0]
	permut=1000
	while float(FETres[0])<0.01 and permut<=100000000:
		FETres = r['fisher.test'](popmatrix,conf_level = 0.95, simulate_p_value=True, B = permut)[0]
		permut=permut*10
	return float(FETres[0])

def FET(x):
	r=robjects.r
	popcounts = robjects.IntVector(x)
	popmatrix = robjects.r['matrix'](popcounts, nrow = 2,byrow=False)
	FETres = r['fisher.test'](popmatrix,conf_level = 0.95)[0]
	return str(FETres[0])
	             
def biallelic(inp,cp):
	h,h1=collections.defaultdict(lambda:0),[]
	for i in cp:
		h["A"]+=int(inp.split()[int(i)+2].split(":")[0])
		h["T"]+=int(inp.split()[int(i)+2].split(":")[1])
		h["C"]+=int(inp.split()[int(i)+2].split(":")[2])
		h["G"]+=int(inp.split()[int(i)+2].split(":")[3])
	while len(h1)<=1:
		h1.append(max(h,key=lambda a: h.get(a)))
		del h[max(h,key=lambda a: h.get(a))]
	return h1
		
def pval_test(x):
	if float(x)<2.2e-16:
		return 2.2e-16
	elif float(x)>=0.9999999999999999:
		return 0.99999999
	else:
		return x

test=[10,100,100,10]
#print FET(test,"less")		

pops=str(options.pops).split(",")

r=robjects.r
stl,gl,ll,c=[],[],[],0
out=open(options.out+".FET","w")
for l in open(options.inp):
	c+=1
	if SNP_test(l,pops,options.maxcov,options.mincov,options.count,options.th)!=0:
		biallele=biallelic(l,pops)
		FETlist,avcov,=[],[]
		for j in pops: 
			FETlist.append(counthash(l,j)[biallele[0]])
			FETlist.append(counthash(l,j)[biallele[1]])				
		print FET(FETlist)
		out.write(l.rstrip()+"\t"+str(FET(FETlist))+"\n")			
	if c%100000==0:
		print str(c)+" positions processed"