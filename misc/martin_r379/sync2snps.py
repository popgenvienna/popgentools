import sys
import collections
import math
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python sync2snps.py -s SNPS.sync > SNPS.snps
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
parser.add_option_group(group)
(options, args) = parser.parse_args()


sync=open(str(options.sync),"r")

#2L	4910	A	6:0:0:0:0:0	-	5:0:0:0:0:0	25:0:0:0:0:0	33:0:0:0:0:0
#2L	4911	G	0:0:0:7:0:0	-	0:0:0:5:0:0	0:0:0:27:0:0	0:0:0:33:0:0
#2L	4912	A	7:0:0:0:0:0	-	5:0:0:0:0:0	27:0:0:0:0:0	32:0:0:0:0:0
#2L	4913	G	0:0:0:7:0:0	-	0:0:0:5:0:0	0:0:0:27:0:0	0:0:0:34:0:0
#2L	4914	A	7:0:0:0:0:0	-	5:0:0:0:0:0	26:0:0:0:0:0	33:0:0:0:0:0
#2L	4915	G	0:0:0:8:0:0	-	0:0:0:5:0:0	0:0:0:22:0:0	0:0:0:33:0:0
#2L	4916	A	8:0:0:0:0:0	-	6:0:0:0:0:0	25:0:0:0:0:0	33:0:0:0:0:0
#2L	4917	G	0:0:0:8:0:0	-	0:0:0:6:0:0	0:0:0:23:0:0	0:0:0:32:0:0
#2L	4918	C	0:0:8:0:0:0	-	0:0:6:0:0:0	0:0:22:0:0:0	0:0:30:0:0:0

print "#chr\tpos\trc\tallele_states\t+"""
dict0={}
dict1={}
for l in sync:
	if l.rstrip("\n")!="":
		a=l.split()
		b=a[0]+"_"+a[1]
		for i in range(0,len(l.split())-3):
			#print i+3
			bases=["A","T","C","G"]
			comp=a[i+3]
			if comp!="-":
				dict0[str(i+3)]=dict(zip(bases,map(float,comp.split(":")[:4])))
		#print dict0
		if a[2]!="A":
			dict1["A"]=sum(value.get("A", 0) for value in dict0.values())
		if a[2]!="T":
			dict1["T"]=sum(value.get("T", 0) for value in dict0.values())
		if a[2]!="C":
			dict1["C"]=sum(value.get("C", 0) for value in dict0.values())
		if a[2]!="G":
			dict1["G"]=sum(value.get("G", 0) for value in dict0.values())
		print a[0]+"\t"+a[1]+"\t"+a[2]+"\t"+max(dict1, key = lambda x: dict1.get(x))+"\t+"
		dict0={}
		dict1={}



