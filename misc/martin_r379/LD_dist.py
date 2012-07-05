import sys
import collections
from optparse import OptionParser,OptionGroup
import copy
import numpy
from statlib.stats import *
from statlib import pstat
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python LD_dist.py -i candidates.cmh -s full.cmh -b 10 -m 1000 -o out.ld
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script takes the set of SNPs in the proximity of candidates SNPs and spits out the average p-values for SNPs in bins until a maximum distance. E.g. -m 500 and -b 20 will result in bins of 50bp length from -500 to +500bp around the candidate.
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-s", "--snps", dest="snps", help="cmh output with all SNPs")
parser.add_option("-b", "--bins", dest="bins", help="number of bins")
parser.add_option("-m", "--maxdist", dest="maxdist", help="maximum distance to candidate.")
parser.add_option("-o", "--out", dest="out", help="outputfile for boxplot")
parser.add_option_group(group)
(options, args) = parser.parse_args()


def posh(x):
	h=collections.defaultdict(lambda:collections.defaultdict(lambda:0.0))
	for l in open(x,"r"):
		if l.rstrip() !="":
			a=l.split()
			chrom=a[0]
			pos=a[1]
			pval=a[-1]
			h[chrom][pos]=pval
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
b=[]
cand=posh(str(options.inp))
full=posh(str(options.snps))
bin=int(options.bins)
maxdist=int(options.maxdist)
min1=-maxdist
for i in range(1,int(bin)+1):
	max1=min1+(maxdist/int(bin)*2)
	b.append(str(min1))
	if min1==-(maxdist/int(bin)*2):
		b.append("cand")
	min1=max1	
	#print b
counter=0
out.write("Chrom\tPos\t"+"\t".join(b)+"\t"+"\t".join(b)+"\n")
for k,v in cand.items():
	for k2,v2 in v.items():
		counter+=1
		count,pav,psd,plist,plit=[],[],[],[],[]
		min1,minc,c=-maxdist,int(k2)-maxdist,0
		for i in range(1,int(bin)+1):
			max1=min1+(maxdist/float(bin)*2)
			maxc=minc+(maxdist/float(bin)*2)
			for k1,v1 in full[k].items():
				if int(k1)>=minc and int(k1)<maxc and int(k1)!=int(k2):
					c+=1
					plist.append(-numpy.log10(float(v1)))
					#plit.append(k1)
			#print k,k2,min1,plit
			#b.append(int(min1))
			count.append(str(c))
			if len(plist)>1:
				pav.append(str(lmean(plist)))
				psd.append(str(lstdev(plist)))
			else:
				pav.append("na")
				psd.append("na")
			c=0
			plist,plit=[],[]
			if min1==-(maxdist/int(bin)*2):
				pav.append(str(-numpy.log10(float(v2))))
				psd.append("0")
			min1=max1
			minc=maxc
		#line(pav,b,str(options.inp))
		out.write(str(k)+"\t"+str(k2)+"\t"+"\t".join(pav)+"\t"+"\t".join(psd)+"\n")
		if counter%2==0:
			print str(counter)+" SNPs processed"