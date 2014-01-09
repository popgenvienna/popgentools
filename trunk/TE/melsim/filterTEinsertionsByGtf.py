import sys
import collections
from optparse import OptionParser, OptionGroup

def loadgtfarray(file):
	"""
	2R      nucmer  orthologous     1300850 1437269 .       .       .       2R:121481-223453 count=29
	2R      nucmer  orthologous     1491008 1499001 .       .       .       2R:226351-234441 count=4
	2R      nucmer  orthologous     1506368 1506648 .       .       .       X:6942-7222 count=1
	2R      nucmer  orthologous     1529046 2148762 .       .       .       2R:237441-840458 count=216
	"""
	gtfa=collections.defaultdict(lambda:[])
	for l in open(file):
		l=l.rstrip("\n")
		a=l.split("\t")
		chr=a[0]
		start=int(a[3])
		end=int(a[4])
		if start > end:
			start,end=(end,start)
		gtfa[chr].append((start,end))
	for chr in gtfa.keys():
		val=gtfa[chr]
		sval=sorted(val,key=lambda v:-(v[1]-v[0]))
	return gtfa

def is_present(chr,pos,gtfa):
	chrspec=gtfa[chr]
	for a in chrspec:
		if pos >= a[0] and pos <=a[1]:
			return True
	return False
	pass

parser = OptionParser()
parser.add_option("--te",dest="te",help="A file containing the TEs to filter")
parser.add_option("--gtf",dest="gtf",help="A file with the gtfs")
(options, args) = parser.parse_args()

gtfa=loadgtfarray(options.gtf)


for l in open(options.te):
	"""
	3R	12259.5	FR	INE-1	0.990340136054422	TIR	-	ncorrdist=177	11933	12032	0.994285714285714	175	174	1	0	12487	12586	0.986394557823129	147	145	2	0
	3R	13787.5	FR	1360	0.936607142857143	TIR	-	ncorrdist=190	13523	13622	0.9375	80	75	5	0	13953	14052	0.935714285714286	14131	9	0
	3R	15517	F	INE-1	0.470149253731343	TIR	-	-	15318	15417	0.470149253731343	134	63	71	0	-	-	-	-	-	-
	3R	15674	R	INE-1	0.556521739130435	TIR	-	-	-	-	-	-	-	-	-	15774	15873	0.556521739130435	115	6451	0
	"""
	l=l.rstrip("\n")
	a=l.split("\t")
	chr=a[0]
	pos=int(float(a[1]))
	ispresent=is_present(chr,pos,gtfa)
	if(ispresent):
		print l