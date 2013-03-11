import sys
import collections
import math
import re
import random
from optparse import OptionParser, OptionGroup

def load_recombination_rate(recfile,minrec):
	recrate={}
	activechr=None
	started=False
	recstart=0
	for l in open(recfile):
		"""
		2L:0..100000            0.00            0.00            0.00      
		2L:100000..200000       0.00            0.00            0.00      
		2L:200000..300000       0.00            0.00            1.89      
		2L:300000..400000       1.89            1.92            1.95      
		2L:400000..500000       1.95            1.98            2.01      
		2L:500000..600000       2.01            2.04            2.07     
		"""
		l=l.rstrip()
		a=l.split("\t")
		m=re.search(r"(\w+):(\d+)..(\d+)",a[0])
		chr=m.group(1)
		start=int(m.group(2))
		end=int(m.group(3))
		totest=((float(a[1]),start),(float(a[2]),(start+end)/2),(float(a[3]),end))
		for t in totest:
			ar=t[0]
			if(not started and ar>=minrec):
				started=True
				recstart=t[1]
			elif(started and ar<minrec):
				started=False
				recrate[chr]=(recstart,t[1])
	return recrate
				
			




parser = OptionParser()
parser.add_option("--input",dest="input",help="A file containing the dgrp haplotypes (MimicrEE input)")
parser.add_option("--recfile",dest="recfile",help="A file containing the low recombining region")
parser.add_option("--recthres",dest="thres",help="Recombination rate threshold")
(options, args) = parser.parse_args()

rr=load_recombination_rate(options.recfile,float(options.thres))

for line in open(options.input):
	line=line.rstrip()
	a=line.split("\t")
	chr=a[0]
	pos=int(a[1])
	rec=rr[chr]
	if(pos >=rec[0] and pos<=rec[1]):
		print line
