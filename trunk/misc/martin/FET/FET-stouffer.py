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
from decimal import *


#Author: Martin Kapun
#version: 2.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python FET_stouffer.py -i candidates.sync -b 100 -p 1,5,9+2,6,10+3,4,8 -s 1,2 -b 100 -c 10 -m 500 -n 10 -o output
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/ 
	""") 


parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-p", "--pops", dest="pops", help="define populations for which you want to calculate the allele frequency change. replicate populations are separated by a \",\" and different treatments or timeposints are separated by a \"+\" (e.g. 1,5,9+2,6,10+3,4,8)")
parser.add_option("-s", "--se", dest="se", help="define start and end of pop comparison of CMH: e.g.: 1,2 for Base-F15, or 1,3 for base-f37")
parser.add_option("-b", "--bins", dest="bins", help="number of bins in histogram")
parser.add_option("-c", "--count", dest="count", help="minor allele count threshold summed across all populations")
parser.add_option("-m", "--maxcov", dest="maxcov", help="maximum coverage threshold for each population separatehd")
parser.add_option("-n", "--mincov", dest="mincov", help="minimum coverage")
parser.add_option("-o", "--out", dest="out", help="output file")


parser.add_option_group(group)
(options, args) = parser.parse_args()

def cov(inp,i): ##returns coverage of population i
	h=collections.defaultdict(lambda:0)
	h["A"]=int(inp.split()[int(i)+2].split(":")[0])
	h["T"]=int(inp.split()[int(i)+2].split(":")[1])
	h["C"]=int(inp.split()[int(i)+2].split(":")[2])
	h["G"]=int(inp.split()[int(i)+2].split(":")[3])
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


def counthash(inp,i): ### returns allelefrequencies fo a given population i
	h=collections.defaultdict(lambda:0)
	h["A"]=int(inp.split()[int(i)+2].split(":")[0])
	h["T"]=int(inp.split()[int(i)+2].split(":")[1])
	h["C"]=int(inp.split()[int(i)+2].split(":")[2])
	h["G"]=int(inp.split()[int(i)+2].split(":")[3])
	return h
	
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

def FET(x,side):
	r=robjects.r
	popcounts = robjects.IntVector(x)
	popmatrix = robjects.r['matrix'](popcounts, nrow = 2,byrow=False)
	FETres = r['fisher.test'](popmatrix,conf_level = 0.95,alternative=side)[0]
	if float(FETres[0])>=0.999999999999:
		return 0.999999999999
	else:
		return float(FETres[0])
	             
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
		
#test=[10,100,100,10]
#print FET(test,"less")
#print FET(test,"less")*1e+56
as_numeric=robjects.r['as.numeric']
getcontext().prec = 50

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
r=robjects.r
stl,gl,ll,c=[],[],[],0
out=open(options.out+".stouffer","w")
for l in open(options.inp):
	c+=1
	if SNP_test(l,allpops,options.maxcov,options.mincov,options.count)!=0:
		fullavcov=0
		biallele=biallelic(l,allpops)
		#print biallele
		pvalhash={"greater":[],"less":[]}
		stoufferhash=collections.defaultdict(lambda:0)
		stouffer=0.0
		for i in range(0,len(comptime)):
			fullavcov+=meanstdv([cov(l,comptime[i][0]),cov(l,comptime[i][1])])[0]
		for i in range(0,len(comptime)):
			FETlist,avcov,=[],[]
			for j in range(0,len(comptime[i])): 
				FETlist.append(counthash(l,comptime[i][j])[biallele[0]])
				FETlist.append(counthash(l,comptime[i][j])[biallele[1]])
				avcov.append(cov(l,comptime[i][j]))	
			for j in sign:
				pvalhash[j].append([FET(FETlist,j),math.sqrt((meanstdv(avcov)[0]/float(fullavcov)))])			
		for tail, replicates in pvalhash.items():
			result=0.0
			for pval,wi in replicates:
				r['options'](digits=22)
				result+=wi*r['qnorm'](pval)[0]
			r['options'](digits=22)
			stouffer=float(r['pnorm'](result)[0])
			if stouffer>=0.999999999999:
				stouffer=0.999999999999
			stoufferhash[tail]=stouffer
			if tail=="greater":
				gl.append(stouffer)
			else:
				ll.append(stouffer)
		minst=min(stoufferhash,key=lambda a: stoufferhash.get(a))	
		out.write(l.rstrip()+"\t"+str(stoufferhash["greater"])+"\t"+str(stoufferhash["less"])+"\t"+minst+"\t"+str(stoufferhash[minst])+"\n")			
		stl.append(stoufferhash[minst])
	if c%10000==0:
		print str(c)+" positions processed"
hist(stl,100,str(options.out)+"_minor_Stouffer","stouffer minor p-value distribution")
hist(gl,100,str(options.out)+"_greater_Stouffer","stouffer p-value distribution of \n one-sided FET, with Ha: greater")
hist(ll,100,str(options.out)+"_less_Stouffer","stouffer p-value distribution of \n one-sided FET, with Ha: less")