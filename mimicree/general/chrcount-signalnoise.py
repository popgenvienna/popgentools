#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections



		



        


def get_selected_hash(file):
	s=set([])
	for line in open(file):
		"""
		2L	6009099	A	0.1	0.5
		3R	12520973	A	0.1	0.5
		"""
		line=line.rstrip()
		a=line.split("\t")
		chr,pos=(a[0],a[1])
		key=chr+":"+pos
		s.add(key)
	return s	

parser = OptionParser()
selected=sys.argv[1]
selhash=get_selected_hash(selected)
filename=sys.argv[2]


selectedsignal=0.0
selectedcount=0.0
normalsignal=0.0
normalcount=0.0
for line in open(filename):
	line=line.rstrip()
	a=line.split("\t")
	chr,pos=(a[0],a[1])
	pval=float(a[-1])
	key=chr+":"+pos
	if key in selhash:
		selectedsignal+=pval
		selectedcount+=1.0
	else:
		normalsignal+=pval
		normalcount+=1.0

signal=selectedsignal/selectedcount
noise=normalsignal/normalcount
sig2noise=signal/noise
print signal,noise,sig2noise
