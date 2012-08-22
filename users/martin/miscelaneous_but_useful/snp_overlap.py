import sys
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="\npython %prog --input1 candidates.fet --input2 input.sync > output.sync"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,

"""
H E L P:
____________

This script uses the first two columns (Chromosome and position) to match both input files (--input1 and --input2) and prints the lines of input2 overlapping with input1
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input1", dest="i", help="any file with Chromosme and positions in the first two columns")
parser.add_option("--input2", dest="j", help="any file with Chromosme and positions in the first two columns")


genehash={}
for l in open(options.i,"r"):
	if len(l.split())>1:
		genehash[l.split()[0]+"_"+l.split()[1]]=l

for l in open(options.j,"r"):
	if len(l.split())>1:
		if l.split()[0]+"_"+l.split()[1] in genehash:
			print l.rstrip()
