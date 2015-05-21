import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import truncsamlib
import math
import re


		



parser = OptionParser()
parser.add_option("--sam",dest="sam",help="the input file as sam")
parser.add_option("--bound",dest="bound",help="the coverage bound")
(options, args) = parser.parse_args()
samh=load_samhash(options.sam)

s=truncsamlib.PTruncSamEntry(options.sam)
avc=s.get_averagecoverage(int(optionbs.bound)
truncs=s.get_truncations()

tc=collections.defaultdict(lambda:0)
for t in truncs:
	tc[t]+=1

print("Av.cov: {0}".format(avc))
for t,c in tc.items():
	freq=float(c)/avc
	print "{0}\t{1}\t{2}\t{3}".format(t[0],t[1],c,freq)

