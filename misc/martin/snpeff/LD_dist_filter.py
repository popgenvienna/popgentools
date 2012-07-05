import sys
import collections
from optparse import OptionParser,OptionGroup
import copy
import numpy
import math
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

1)	usage: python LD_dist.py -i candidates.cmh -s full.cmh -p 200 -m 2 -t 0.05 -o out.ld
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script takes the set of SNPs in the proximity of candidates SNPs and spits out the average p-values for SNPs in bins until a maximum distance. E.g. -m 500 and -b 20 will result in bins of 50bp length from -500 to +500bp around the candidate.
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-s", "--snps", dest="snps", help="cmh output with all SNPs")
parser.add_option("-p", "--step", dest="step", help="number of bins")
parser.add_option("-m", "--minsnps", dest="minsnps", help="minimum number of SNPs in bins")
parser.add_option("-o", "--out", dest="out", help="outputfile for boxplot")
parser.add_option("-t", "--th", dest="th", help="threshold")
parser.add_option_group(group)
(options, args) = parser.parse_args()

def meanstdv(x): ### calculate mean, median, stdev. standard error: x=datalist
    n,mean,std,se,median = len(x),0,0,0,0
    for a in x:
        mean = mean + a
    mean = mean / float(n)
    return mean

def posh(x):
    h=collections.defaultdict(lambda:collections.defaultdict(lambda:0.0))
    for l in open(x,"r"):
        if l.rstrip() !="":
            a=l.split()
            chrom=a[0]
            pos=int(a[1])
            pval=float(a[-1])
            h[chrom][pos]=-math.log10(pval)
    return h

def nbh(x):
    has=collections.defaultdict(lambda:collections.defaultdict(lambda:[]))
    a=0
    b=0
    for key,iz in x.items():
        for item in sorted(map(int,iz.keys())):
            has[key][b]=[a,item]
            a=copy.deepcopy(b)
            b=copy.deepcopy(item)
        has[key][b]=[a,0]
        a=0
        b=0
    return has

out=open(str(options.out),"w") #### select output file 
cand=posh(str(options.inp))   ###produce hash of cand snps in the format {chrom:pos:-log10(pval)}
candout=posh(str(options.inp)) ###produce hash of cand in the format {chrom:pos{-log10(pval)}}
full=posh(str(options.snps)) ###produce hash of all snps in the format {chrom:pos{-log10(pval)}}
print "loading full SNP datatset done"
fullsnps=nbh(full) ### produce a hash in the format {chrom:pos:[up_SNP,down_SNP]}
print "calculating adjacencies done"
step=int(options.step)
minsnps=int(options.minsnps)
min1=0
loci=collections.defaultdict(lambda:[])
th=float(options.th)
counter=0
for l in open(str(options.inp),"r"):         #### loop through all chromosomal arms
    a=l.split()   ### loop through all positions on chromsomal arms
    k=a[0]
    k2=int(a[1])
    v2=-math.log10(float(a[-1]))
    #print l.rstrip()
    if k2 in candout[k]:
        counter+=1
        bin=collections.defaultdict(lambda:0.0)
        bin[k2]=v2
        bin1=copy.deepcopy(bin)
        lim=fullsnps[k][min(bin.keys())][0]
        count=1
        while count<minsnps and lim in full[k] and lim!=0:
            count+=1
            lim=fullsnps[k][min(bin.keys())][0]
            bin[lim]=full[k][lim]
        newp=float(meanstdv(bin.values()))
        pth=float(v2)/float(minsnps)
        length=max(bin.keys())-min(bin.keys())
        #print "up",bin,k2,lim,pth,newp,length
        while newp>=pth and length<=step and pth>-math.log10(th):
            bin1=copy.deepcopy(bin)
            lim=fullsnps[k][min(bin.keys())][0]
            if lim!=0:
                bin[lim]=full[k][lim]
                minsnps+=1
                newp=float(meanstdv(bin.values()))
                pth=float(v2)/float(minsnps)
                #print "up",bin,k2,lim,pth,newp,length
            else:
                break
        loci[l].append(k2-min(bin1.keys()))
        loci[l].append(minsnps)
        for can,vol in bin1.items():
            if can in candout[k]:
                print "del",k,can
                del candout[k][can]
                
                
        minsnps=int(options.minsnps)
        bin=collections.defaultdict(lambda:0.0)
        bin[k2]=v2
        bin1=copy.deepcopy(bin)
        lim=fullsnps[k][max(bin.keys())][1]
        count=1
        while count<minsnps and lim in full[k] and lim!=0:
            count+=1
            lim=fullsnps[k][max(bin.keys())][1]
            bin[lim]=full[k][lim]
        newp=float(meanstdv(bin.values()))
        pth=float(v2)/float(minsnps)
        length=max(bin.keys())-min(bin.keys())
        #print "down",bin,k2,lim,pth,newp,length
        while newp>=pth and length<=step and pth>-math.log10(th):
            bin1=copy.deepcopy(bin)
            lim=fullsnps[k][max(bin.keys())][1]
            if lim!=0:
                bin[lim]=full[k][lim]
                minsnps+=1
                newp=float(meanstdv(bin.values()))
                pth=float(v2)/float(minsnps)
                #print "down",bin,k2,lim,pth,newp,length
            else:
                break
        loci[l].append(max(bin1.keys())-k2)
        loci[l].append(minsnps)
        for can,vol in bin1.items():
            if can in candout[k]:
                print "del",k,can
                del candout[k][can]
                
        if counter%10==0:
            print str(counter)+" SNPs processed"
for k,v in loci.items():
    out.write(k.rstrip()+"\t"+"\t".join(map(str,v))+"\t"+str(+v[0]+v[2])+"\n")