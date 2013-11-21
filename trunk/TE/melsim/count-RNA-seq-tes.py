import sys
import random
from optparse import OptionParser, OptionGroup
import collections




def getTEtranslator(file,targethierlevel):
	"""
	KEPLER	FBgn0063570	Dbuz\Kepler	Foldback	Foldback	DNA
	1360	FBgn0005673	1360	SINE	non-LTR	RNA
	DNTOMRETA	FBgn0004357	Dana\Tom	LTR	LTR	RNA
	PPI251	FBgn0003055	P-element	SINE	non-LTR	RNA
	"""
	hier={}
	targetindex=0
	for line in open(file):
		line=line.rstrip("\n")
		a=line.split("\t")

		if(line.startswith("insert")):
			found=False
			for i,e in enumerate(a):
				if e==targethierlevel:
					targetindex=i
					found=True
					break
			if not found:
				raise ValueError("Could not find hierarchy level "+targethierlevel + " in header: "+line)
		else:
			totranslate=a[0]
			intotranslate=a[targetindex]
			if totranslate in hier:
				raise ValueError("Invalid; " +id+" present multiple times")
			hier[totranslate]=intotranslate
	return hier

parser = OptionParser()
parser.add_option("--input",		dest="sam",help="A sam file")
parser.add_option("--hierarchy",	dest="hier",help="A hierarchy of TE sequences")
parser.add_option("--hierarchy-level",	dest="level",help="Level of the targeted TE hierarchy")
parser.add_option("--min-map-qual",	dest="minmapqual",help="Minimum mapping quality")
(options, args) = parser.parse_args()
mmq=int(options.minmapqual)

translator=getTEtranslator(options.hier,options.level)
mapcount=0
maptecount=0
maptewithqual=0
testat=collections.defaultdict(lambda:0)
for line in open(options.sam):
	# FCD0KN9ACXX:8:1101:1228:2236#NNNNNNNN   16      3R      27193650        170     100M    *       0       0       ATAGGCTGNNATCGATGCTGGACTTTCAATGCTGCAGTGTGTGCGTCAGCTGGACAGACAAGAAAAACAATTAGGCAGTAGTACCAGCCCTAAAACCCTC    bccc_WKKBBcccccbccccdbddddeeeeeeggggghihhhiihhhiiiiihhiiigiiiiiiihfiiiiiiiiiiiiihgeiiiigggggecceea_a    AS:i:92 XS:i:0  XF:i:3  XE:i:2  NM:i:1
	line=line.rstrip("\n")
	if line.startswith("@"):
		continue
	a=line.split("\t")
	flag=int(a[1])
	if flag& 0x004:
		continue # discard unmapped
	mapcount+=1
	chr=a[2]
	mq=int(a[4])
	if chr in translator:
		maptecount+=1
		if mq < mmq:
			continue
		maptewithqual+=1
		fam=translator[chr]
		testat[fam]+=1
it=sorted(testat.items(),key=lambda k:-k[1])


print "Reads mapped\t{0}".format(mapcount)
print "Reads mapped to TE\t{0}".format(maptecount)
print "Reads mapped to TE with quality of {0}\t{1}".format(mmq,maptewithqual)
for t in it:
	print "{0}\t{1}".format(t[0],t[1])

	