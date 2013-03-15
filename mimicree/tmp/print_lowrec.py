#!/usr/bin/env python

import sys
import collections
import math
import re
import random
from optparse import OptionParser, OptionGroup

class RecombinationWindow:
	def __init__(self,chr,start,end,rec):
		self.chr=chr
		self.start=start
		self.end=end
		self.rec=rec




def load_recombination(recfile,minrec):
	recwindows=load_recombinaton_file(recfile)
	recchr=group_by_chromosome(recwindows)
	
	recrate={}
	for chr,windows in recchr.items():
		started=False
		recstart=0
		for w in windows:
			if(not started and w.rec >= minrec):
				started=True
				recstart=w.start
			elif(started and w.rec<minrec):
				recrate[w.chr]=(recstart,w.start)
				break
	assert(len(recrate.keys())==5) # 5 chromosomes of dmel
	return recrate
		
	

def group_by_chromosome(recwindows):
	toret=collections.defaultdict(lambda:[])
	for w in recwindows:
		toret[w.chr].append(w)
	return toret

def load_recombinaton_file(recfile):
	toret=[]
	for l in open(recfile):
		"""
		2L:0..100000            0.00            0.00            0.00      
		2L:100000..200000       0.00            0.00            0.00      
		2L:200000..300000       0.00            0.00            1.89      
		2L:300000..400000       1.89            1.92            1.95      
		2L:400000..500000       1.95            1.98            2.01      
		2L:500000..600000       2.01            2.04            2.07
		2L:600000..700000       1.0             0.00            0.00     
		"""
		l=l.rstrip()
		a=l.split("\t")
		m=re.search(r"(\w+):(\d+)..(\d+)",a[0])
		chr=m.group(1)
		start=int(m.group(2))
		end=int(m.group(3))
		recrate=float(a[2])
		
		rw=RecombinationWindow(chr,start,end,recrate)
		toret.append(rw)
	return toret

parser = OptionParser()
parser.add_option("--recfile",dest="recfile",help="A file containing the low recombining region")
parser.add_option("--recthres",dest="thres",help="Recombination rate threshold")
(options, args) = parser.parse_args()

rr=load_recombination(options.recfile,float(options.thres))
print rr
