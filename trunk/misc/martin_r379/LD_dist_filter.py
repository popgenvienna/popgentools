import sys
import collections
from optparse import OptionParser,OptionGroup
import copy
import numpy
#from statlib.stats import *
#from statlib import pstat
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python LD_dist.py -i candidates.cmh -s full.cmh -p 200 -o out.ld
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script takes the set of SNPs in the proximity of candidates SNPs and spits out the average p-values for SNPs in bins until a maximum distance. E.g. -m 500 and -b 20 will result in bins of 50bp length from -500 to +500bp around the candidate.
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-s", "--snps", dest="snps", help="cmh output with all SNPs")
parser.add_option("-p", "--step", dest="step", help="number of bins")
parser.add_option("-m", "--minsnps", dest="minsnps", help="minimum number of SNPs in bins")
parser.add_option("-o", "--out", dest="out", help="outputfile for boxplot")
parser.add_option_group(group)
(options, args) = parser.parse_args()

def meanstdv(x): ### calculate mean, median, stdev. standard error : x=datalist
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

def posh(x):
	h=collections.defaultdict(lambda:collections.defaultdict(lambda:0.0))
	for l in open(x,"r"):
		if l.rstrip() !="":
			a=l.split()
			chrom=a[0]
			pos=a[1]
			pval=a[-1]
			h[chrom][str(pos)]=float(pval)
	return h


def sort_dict(dictionary,choice,ct):
	h={}
	k=dictionary.keys()
	keys=map(int,k)
	if choice=="reverse":
		be=sorted(keys,reverse=True)
		ct1=0
		for a in be:
			if ct1<ct: 
				h[a]=dictionary[a]
			ct1+=1
	else:
		be=sorted(keys)
		ct1=0
		for a in be:
			if ct1<ct: 
				h[a]=dictionary[a]
			ct1+=1
	return h
	

def line(x,bins,out):
	import rpy2.robjects as robjects
	from rpy2.robjects.packages import importr
	r=robjects.r
	grdevices = importr('grDevices')
	#rang=robjects.vectors.FloatVector([-4.75,-4.74])
	values=robjects.vectors.FloatVector(x)
	xaxis=robjects.vectors.FloatVector(bins)
	xlim=robjects.vectors.FloatVector([-1100,1100])
	grdevices.png(out+'_line.png')
	r.plot(xaxis,values,ylab="log[p-value]",xlab="bins",type="o",)
	r.axis(1,tick=-0.001)
	grdevices.dev_off()
	
out=open(str(options.out),"w")
cand=posh(str(options.inp))
candout=posh(str(options.inp))
full=posh(str(options.snps))
disthash=collections.defaultdict(lambda:collections.defaultdict(lambda:[]))
counter=0
#print candout
step=int(options.step)
minsnps=int(options.minsnps)
min1=0
for l in open(str(options.inp),"r"):         #### loop through all chromosomal arms
	counter+=1
	a=l.split()   ### loop through all positions on chromsomal arms
	k=a[0]
	k2=a[1]
	v2=a[-1]
	if str(a[1]) in candout[k]:
		contind,continu=1,1
		fullbin,bin2,bl,length=collections.defaultdict(lambda:0),[],0,0                       
		for k1,v1 in full[k].items():  ### loop through all SNPs
			if int(k1)>=int(k2)-step and int(k1)<int(k2):  ### consider if position smaller than candidate and larger than candidate - Stepsize
				fullbin[int(k1)]=-numpy.log10(float(v1)) #### key: position, value -log10(pvalue)
		if len(fullbin)>=minsnps:  ### number of entries in dictionary larger or equal to limit??
			if float(meanstdv(sort_dict(fullbin,"reverse",minsnps).values())[0])>=float(-numpy.log10(float(v2))/2):
				bin=sort_dict(fullbin,"reverse",minsnps)
				for k3,v3 in bin.items():
					if str(k3) in candout[k]:
						#print k3
						del candout[k][str(k3)]
				lim=min(sort_dict(bin,"reverse",minsnps).keys())
				while continu==1:  ### continue if limits not reached 

					for k4 in full[k].keys():
						if int(k4)<int(lim): 
							bin2.append(int(k4))
					bin1=sorted(bin2)
					bin[bin1[-1]]=float(-numpy.log10(full[k][str(bin1[-1])]))
					min1=copy.deepcopy(minsnps)
					minsnps+=1
					if float(meanstdv(sort_dict(bin,"reverse",minsnps).values())[0])>=float(-numpy.log10(float(v2))/2):
						if str(bin1[-1]) in candout[k]:
							print "delete", bin1[-1]
							del candout[k][str(bin1[-1])]
						bin2=[]
						lim=min(bin.keys())
					else:
						disthash[k][a[1]].append(str(int(a[1])-int(lim)))
						disthash[k][a[1]].append(str(min1))
						bin2=[]
						continu=0
						minsnps=int(options.minsnps)
						
			else:
				disthash[k][a[1]].append("0")
				disthash[k][a[1]].append(str(minsnps))
				minsnps=int(options.minsnps)
		else:
			disthash[k][a[1]].append("0")
			disthash[k][a[1]].append("0")
			minsnps=int(options.minsnps)

		fullbin,bin2,bl,length=collections.defaultdict(lambda:0),[],0,0
		for k1,v1 in full[k].items():  ### loop through all SNPs
			if int(k1)<=int(k2)+step and int(k1)>int(k2):  ### consider if position smaller than candidate and larger than candidate - Stepsize
				fullbin[int(k1)]=-numpy.log10(float(v1)) #### key: position, value -log10(pvalue)
		if len(fullbin)>=minsnps:  ### number of entries in dictionary larger or equal to limit??
			if float(meanstdv(sort_dict(fullbin,"",minsnps).values())[0])>=float(-numpy.log10(float(v2))/2):
				bin=sort_dict(fullbin,"",minsnps)
				for k3,v3 in bin.items():
					if str(k3) in candout[k]:
						print "delete", k3
						del candout[k][str(k3)]
				lim=max(sort_dict(bin,"",minsnps).keys())
				while contind==1:      ### continue if limits not reached       
					for k4 in full[k].keys():
						if int(k4)>int(lim): 
							bin2.append(int(k4))
					bin1=sorted(bin2)

					bin[bin1[0]]=float(-numpy.log10(full[k][str(bin1[0])]))
					min1=copy.deepcopy(minsnps)
					minsnps+=1     

					if float(meanstdv(sort_dict(bin,"",minsnps).values())[0])>=float(-numpy.log10(float(v2))/2):
						if str(bin1[0]) in candout[k]:
							print "delete", bin1[0]
							del candout[k][str(bin1[0])]
						bin2=[]
						lim=max(bin.keys())
					else:
						disthash[k][a[1]].append(str(int(lim)-int(a[1])))
						disthash[k][a[1]].append(str(min1))
						disthash[k][a[1]].append(str(v2))
						bin2=[]
						contind=0
						minsnps=int(options.minsnps)
						
			else:
				disthash[k][a[1]].append("0")
				disthash[k][a[1]].append("0")
				disthash[k][a[1]].append(str(v2))
				minsnps=int(options.minsnps)
		else:
			disthash[k][a[1]].append("0")
			disthash[k][a[1]].append("0")
			disthash[k][a[1]].append(str(v2))
			minsnps=int(options.minsnps) 
			
	if counter%10==0:
		print str(counter)+" SNPs processed"
for k,v in candout.items():
	for k1,v1 in v.items():
		 out.write(str(k)+"\t"+str(int(k1)-1)+"\t"+str(k1)+"\tfeature\t"+"\t".join(disthash[k][k1][:-1])+"\t"+str(int(disthash[k][k1][0])+int(disthash[k][k1][2]))+"\t"+str(int(disthash[k][k1][1])+int(disthash[k][k1][3]))+"\t"+disthash[k][k1][4]+"\n")	 


