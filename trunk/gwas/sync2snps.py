import sys
import collections
import math
import warnings
from optparse import OptionParser,OptionGroup
from modules.CMH import SyncReader

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python sync2snpeff.py -s SNPS.sync > SNPS.snps
2)	This script prints the reference and the alternative Allele to a list. If the SNP has more than two states it will print the allele with the higher frequency. It takes a sync file (output form synchronize_pileups.pl) as an input which has do be filtered to contain only SNPs (e.g. which filter_sync.py). See output example:

#chr    pos     rc      allele_states   +
2L      5002    G       T       +
2L      5017    T       C       +
2L      5437    A       T       +
2L      5465    C       A       +
2L      5495    A       T       +
2L      5750    A       T       +
2L      5762    T       C       +
2L      5813    G       T       +
2L      5845    C       A       +
2L      6079    C       T       +
2L      6445    T       G       +
2L      6631    A       G       +
2L      6686    A       T       +
2L      6690    C       G       +
2L      6764    T       A       +
2L      6921    G       T       +
2L      6930    A       T       +
2L      6931    C       G       +
2L      7088    A       T       +

	c1:	chromsome
	c2:	SNP position
	c3:	reference allele
	c4:	alternative allele
	c5:	SNP orientation (this is a column required by SNPeff. I do not fully understand how to infer the orinetation of a SNP and thus always put a \"+\" here. Please correct me if I am wrong!!!)
	""")

#########################################################   CODE   #########################################################################

parser.add_option("-s", "--sync", dest="sync", help="sync file of single chromosomal arm")
parser.add_option("--test",action="store_true", dest="test",help="run the doctest")
parser.add_option_group(group)
(options, args) = parser.parse_args()

def maxkeyfornotrefc(dict,refc):
	"""
	>>> maxkeyfornotrefc({"A":5,"T":3,"C":2,"G":1},"T")
	'A'
	>>> maxkeyfornotrefc({"A":5,"T":3,"C":2,"G":1},"A")
	'T'
	>>> maxkeyfornotrefc({"A":5,"T":3,"C":8,"G":1},"A")
	'C'
	>>> maxkeyfornotrefc({"A":5,"T":3,"C":8,"G":1},"C")
	'A'
	>>> maxkeyfornotrefc({"A":5,"T":3,"C":8,"G":1},"N")
	'C'
	
	In case when the second highest are equal what to do then
	>>> maxkeyfornotrefc({"A":5,"T":5,"C":8,"G":1},"C")
	'A'
	
	In case when not actualy a SNP
	>>> maxkeyfornotrefc({"A":5,"T":0,"C":0,"G":0},"A")
	'N'
	"""
	ar=[ (key,value)for key,value in dict.items()]
	ar.sort(key=lambda x: -x[1])
	if(ar[0][0]==refc):
		ar.pop(0)
	secondchar=ar[0][0]
	secondcount=ar[0][1]
	if secondcount==0:
		warnings.warn("Sync entry is no SNP")
		secondchar="N"
	return secondchar


if(options.test):
	import doctest
	doctest.testmod(verbose=1)
	sys.exit()

#2L	4910	A	6:0:0:0:0:0	-	5:0:0:0:0:0	25:0:0:0:0:0	33:0:0:0:0:0
#2L	4911	G	0:0:0:7:0:0	-	0:0:0:5:0:0	0:0:0:27:0:0	0:0:0:33:0:0
#2L	4912	A	7:0:0:0:0:0	-	5:0:0:0:0:0	27:0:0:0:0:0	32:0:0:0:0:0
#2L	4913	G	0:0:0:7:0:0	-	0:0:0:5:0:0	0:0:0:27:0:0	0:0:0:34:0:0
#2L	4914	A	7:0:0:0:0:0	-	5:0:0:0:0:0	26:0:0:0:0:0	33:0:0:0:0:0
#2L	4915	G	0:0:0:8:0:0	-	0:0:0:5:0:0	0:0:0:22:0:0	0:0:0:33:0:0
#2L	4916	A	8:0:0:0:0:0	-	6:0:0:0:0:0	25:0:0:0:0:0	33:0:0:0:0:0
#2L	4917	G	0:0:0:8:0:0	-	0:0:0:6:0:0	0:0:0:23:0:0	0:0:0:32:0:0
#2L	4918	C	0:0:8:0:0:0	-	0:0:6:0:0:0	0:0:22:0:0:0	0:0:30:0:0:0

print "#chr\tpos\trc\tallele_states\t+"
for p in SyncReader(options.sync):
	b=str(p.chr)+"_"+str(p.pos)
	bh={"A":0,"C":0,"G":0,"T":0}
	for pop in p.populations:
		bh["A"]+=pop.A
		bh["T"]+=pop.T
		bh["C"]+=pop.C
		bh["G"]+=pop.G
	maxnonrefc=maxkeyfornotrefc(bh,p.refc)
	print p.chr+"\t"+str(p.pos)+"\t"+p.refc+"\t"+maxnonrefc+"\t+"




