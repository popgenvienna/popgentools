import sys
import collections
from optparse import OptionParser,OptionGroup
import copy
import math

#Author: Martin Kapun & Robert Kofler
#########################################################   HELP   #########################################################################

usage="python %prog --candidates candidates.cmh --snps full.cmh --bin-size 10 --maxdist 1000 --measure median > candidates.ld"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,"""				
H E L P:
_________

Description:
This script uses sync files with any statistic (P-value, allele frequency change, FST,etc.) in the last column as input. Based on the parameters --bin-size and --maxdist, the script averages the statistic of SNPs around all candidates in bins of a given length until a given distance.  Note that the script by default log-scales the statistic, which only applies to P-Values. Thus if you use any other statistic disable log-scaling by setting --log False. Additionally you can choose between calcultating the geometric mean and the median with the parameter --measure. Note that the candidate is excluded from the data. The output consists of two columns: The starting point of the bin relative to the center (the candidate position) and the averaged statistic 
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--candidates", dest="cand", help="*.sync or cmh output file")
parser.add_option("--snps", dest="snps", help="cmh output with all SNPs")
parser.add_option("--bin-size", dest="binsize", help="the size of the bins")
parser.add_option("--distance", dest="maxdist", help="the distance to either side of a candidate SNP")
parser.add_option("--measure", dest="measure", help="What should be calculated, median (median) or geometric mean (gm)")
parser.add_option("--log", dest="log", help="use -log10??, Boolean!!", default=True)


def load_candidates(candidateFile):
	cch=collections.defaultdict(lambda:[])
	for line in open(candidateFile):
		line=line.rstrip()
		a=line.split('\t')
		chr=a[0]
		pos=int(a[1])
		cch[chr].append(pos)
	return cch
		
def median(x):
	"""
	>>> median([float("1e-12"),float("1e-14"),float("1e-10")])
	9.9999999999999998e-13
	>>> median([float("1e-11"),float("1e-20"),float("1e-12")])
	9.9999999999999998e-13
	>>> median([float("1e-8"),float("1e-11"),float("1e-12"),float("1e-20")])
	3.1622776601683786e-12
	>>> median([float("1e-8"),float("1e-10")])
	1.0000000000000007e-09
	"""
	mid =int(len(x)/2)
	sort=sorted(x)
	if len(x)==0:
		return None
	if len(x)%2==0:
		lower=sort[mid-1]
		upper=sort[mid]
		return (float(lower)+float(upper))/2.0
	else:
		return sort[mid]

def mean(x): ### calculate mean, median, stdev. standard error : x=datalist
	from math import sqrt
	n,mean = len(x),0
	for a in x:
		mean += a
	mean = mean / float(n)
	return mean




# median, and log=true
parser.add_option_group(group)
(options, args) = parser.parse_args()
binsize=int(options.binsize)
distance=int(options.maxdist)
bins=math.ceil((2*distance)/binsize)



stat=collections.defaultdict(lambda:[])
# chromosome candidate hash
cch=load_candidates(options.cand)


# iterate of SNPs
for line in open(options.snps):
	line=line.rstrip();
	a=line.split("\t")
	chr=a[0]
	pos=int(a[1])
	if options.log==True:
		pval=-math.log10(float(a[-1]))
	else:
		pval=float(a[-1])
	
	ch=cch[chr]
	for candpos in ch:
		rangestart=candpos-distance
		rangeend=candpos+distance-1
		if(pos<rangestart):
			continue
		if(pos>rangeend):
			continue
		if(candpos==pos): # get rid of the candidate!
			continue
		relativeposition=pos-rangestart
		binindex=int(relativeposition/binsize)
		stat[binindex].append(pval)
start=-distance
end=distance
interval=start

for i in sorted(stat.keys()):
	ar=stat[i]
	if options.measure==median:
		med=median(ar)
	else:
		med=mean(ar)
	print "{0}\t{1}".format(interval,med)
	interval+=binsize


