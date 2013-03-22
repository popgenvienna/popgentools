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


chrguidefile=sys.argv[1]

selectedar=[] # content will be sets with the key: "chr:pos" as string
for i in range(2,len(sys.argv)):
	filename=sys.argv[i]
	selectedar.append(get_selected_hash(filename))
	
for line in open(chrguidefile):
	line=line.rstrip()
	a=line.split("\t")
	chr,pos=(a[0],a[1])
	key=chr+":"+pos
	topr=[]
	for s in selectedar:
		if key in s:
			topr.append("1")
		else:
			topr.append("0")
	toprstr="\t".join(topr)
	print toprstr
	
	
	
	
	