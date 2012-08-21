import sys
import collections
import math
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
usage="python %prog --sync input.sync --snps subset.snps > input.vcf"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'''		

H E L P :
_________

This script takes a file with polymorphic SNPs (e.g. the output of The CMH test) in sync file format as input and identifies the reference and an alternative allele (in case of more than two alleles, the two most common ones will be used) across all populations in the file. Then it creates an output in the VCF v.4 output format, which can be used as the input for SNPeff, etc. Note that only the columns REF and ALT in the VCF file will be filled. all other columns will be left blank (with a "." ). If you only want to convert a subset of the SNP input to vcf you can use the parameter --snps and provide an input file, which has to have at least two columns, specifying the Chromosome and the positions. SNPs in this file need to be a subset of the --sync input. Per default, this option is set to "no" and does not need to be set.

See output example:

##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
3R	4704124	.	T	C	.	.	.
3R	8566126	.	A	C	.	.	.
3R	1178449	.	C	T	.	.	.
3R	15220143	.	C	T	.	.	.
3R	13777597	.	G	A	.	.	.

	''')

#########################################################   CODE   #########################################################################

parser.add_option("-s", "--sync", dest="sync", help="sync file of single chromosomal arm")
parser.add_option("-p", "--snps", dest="snps", help="define a subset of the SNPs which should be converted to VCF",default="no")
parser.add_option_group(group)
(options, args) = parser.parse_args()

snplist=collections.defaultdict(lambda:0)
if options.snps=="no":
	file=options.sync
else: 
	file=options.snps
	
for l in open(file,"r"):	
	a=l.split()
	b=a[0]+"_"+a[1]
	snplist[b]
sync=open(str(options.sync),"r")

print "##fileformat=VCFv4.0"
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
dict0={}
dict1={}
for l in sync:
	if l.rstrip("\n")!="":
		a=l.split()
		b=a[0]+"_"+a[1]
		if b in snplist:
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
			print a[0]+"\t"+a[1]+"\t.\t"+a[2]+"\t"+max(dict1, key = lambda x: dict1.get(x))+"\t.\t.\t."
			dict0={}
			dict1={}



