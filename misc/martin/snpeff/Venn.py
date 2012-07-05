import sys
from optparse import OptionParser, OptionGroup
import collections
import copy


#Author: Martin Kapun
#version 1.0

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python VENN.py -i BF15.igv,BF37.igv,F15F37.igv -n BF15,BF37,F15F37 > snps.venn
2)	script produces presence absence table of all SNPs in the datasets in --i
	""")
	
	
#########################################################   CODE   #########################################################################

parser.add_option("-n", "--names", dest="names", help="names of datasets")
parser.add_option("-i", "--inp", dest="inp", help="input file with candidate SNPs, separated by a \",\"")
parser.add_option_group(group)
(options, args) = parser.parse_args()


igvs=str(options.inp).split(",")
names=str(options.names).split(",")
name={}.fromkeys(str(options.names).split(","),"FALSE")
snps=collections.defaultdict(lambda:copy.deepcopy(name))
for i in range(len(igvs)):
	for l in open(igvs[i],"r"):
		if l.rstrip!="" and "Chromosome" not in l:
			a=l.split()
			snps[a[0]+"_"+a[1]]

for i in range(len(igvs)):
	for l in open(igvs[i],"r"):
		if l.rstrip!="" and "Chromosome" not in l:
			a=l.split()
			snps[a[0]+"_"+a[1]][names[i]]="TRUE"

#print snps
print 	"\t".join(names)
for k,v in snps.items():
	print "\t".join(v.values())