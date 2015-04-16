from gtfIO import GTFReader,GTFWriter;
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import re

parser = OptionParser()
parser.add_option("--input",dest="input",help="A gtf file containing a TE annotation")
(options, args) = parser.parse_args()

def get_fam(comment):
	a=comment.split(" ") 
	target=a[1]
	target=target.lstrip('"Motif:')
	target=target.rstrip('"')
	return target


lengh=collections.defaultdict(lambda:0)
famh=collections.defaultdict(lambda:0)
for e in GTFReader(options.input):
	leng=(e.end-e.start)+1
	lengh[leng]+=1
	fam=get_fam(e.comment)
	famh[fam]+=1

for l in sorted(lengh.keys()):
	print "l\t{0}\t{1}".format(l,lengh[l])

for f in famh.keys():
	print "f\t{0}\t{1}".format(f,famh[f])








