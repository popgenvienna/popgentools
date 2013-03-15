import sys
import collections
import math
import re
import random
from optparse import OptionParser, OptionGroup

	
def get_recrate():
	# manual recrate for dmel. threshold < 1.0
	rr={'X':(1500000,20800000),'2L':(300000,16600000),'2R':(3900000,20700000),'3L':(900000,17400000),'3R':(6600000,25700000)}
	return rr
	# Control, wrote a script with threshold 1.0
	# script output=
	#{'2L': (300000, 16600000), '3R': (6600000, 25700000), 'X': (1500000, 20800000), '2R': (3900000, 20700000), '3L': (900000, 17400000)}
			




parser = OptionParser()
parser.add_option("--input",dest="input",help="A file containing the dgrp haplotypes (MimicrEE input)")
#parser.add_option("--recfile",dest="recfile",help="A file containing the low recombining region")
#parser.add_option("--recthres",dest="thres",help="Recombination rate threshold")
(options, args) = parser.parse_args()

rr=get_recrate();

for line in open(options.input):
	line=line.rstrip()
	a=line.split("\t")
	chr=a[0]
	pos=int(a[1])
	rec=rr[chr]
	if(pos >=rec[0] and pos<=rec[1]):
		print line
