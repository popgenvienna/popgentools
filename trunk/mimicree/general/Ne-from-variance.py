import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import scipy as sp

N=1000
def boost_toobservation(a):
	toret=[]
	for i,count in enumerate(a):
		for k in range(0,int(count)):
			toret.append(i)
	return toret

for line in open(sys.argv[1]):
	"""
	matdist 1       1       164     243     250     198     92      35      12      3       0       3
	matdist 1       2       134     270     269     179     102     33      9       4
	matdist 1       3       129     278     259     190     102     30      9       2       1
	matdist 1       4       143     291     239     165     104     38      14      4       1       0       1
	"""
	if(not line.startswith("matdist")):
		continue
	
	line=line.rstrip()
	a=line.split("\t")
	a.pop(0)
	replicate=a.pop(0)
	generation=a.pop(0)
	obs=boost_toobservation(a)
	mean=sp.mean(obs)
	variance=sp.var(obs)
	
	
	
	ne=(4*N)/(variance+mean**2-mean)
	topr=[]
	topr.append(replicate)
	topr.append(generation)
	topr.append(str(ne))
	toprstr="\t".join(topr)
	print(toprstr)
	