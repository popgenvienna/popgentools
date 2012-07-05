import sys
import re
import collections
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: gff2fasta.py -f genome.fasta -g exonerate.gff --c /PopGenTools/syn-nonsyn/codon-table.txt -o genes 
2)  This script extracts the coding part of the exonerate.gff3 file from the genomic fasta file. It produces two outputs: 1) all genes translated in proteins and 2) all genes as nucleotides. Genes which occur more than once will have a \"_1\" appendeded to the genename to distinguish them. The header of the multifasta output contains the genename and the chromosome as well as the gene coordinates according to the gff file.
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-f", "--fasta", dest="fasta", help="genome")
parser.add_option("-g", "--gff", dest="gff", help="exonerate gff")
parser.add_option("-o", "--out", dest="out", help="outputfile containing SNPs, which show high freq changes starting from very low (0) frequencies")
parser.add_option("-c", "--codon", dest="codon", help="codon table, e.g. /PopGenTools/syn-nonsyn/codon-table.txt")
parser.add_option_group(group)
(options, args) = parser.parse_args()

fastahash={}
genehash,ghash,codon,prothash,infohash=collections.defaultdict(lambda:""),collections.defaultdict(lambda:""),collections.defaultdict(lambda:""),collections.defaultdict(lambda:""),collections.defaultdict(lambda:"")
contigs=[]
def complement(s):
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
	letters = list(s) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters)
def reversecomp(s):
	s=reverse(s)
	s=complement(s)
	return s
for l in open(options.fasta,"r"):
	if ">" in l:
		fastahash[l[1:].split()[0]]=""
		a=l[1:].rstrip()
	if ">" not in l:
		fastahash[a]+=l.rstrip()

for l in open(options.codon,"r"):
	if "#" not in l:
		a=l.split(": ")
		codon[a[0].rstrip()]=a[1].rstrip()

for l in open(options.gff,"r"):
	if l.split()[0] in fastahash:
		if l.split()[2]=="gene":
			reg=re.search(r"ID=(.*);Name.*",l.split()[8])
			if reg.group(1) not in infohash:
				infohash[reg.group(1)]=reg.group(1)+"\t"+str(l.split()[0])+"\t"+str(l.split()[3])+"."+str(l.split()[4])
			if reg.group(1) in infohash:
				infohash[reg.group(1)+"_1"]=reg.group(1)+"\t"+str(l.split()[0])+"\t"+str(l.split()[3])+"."+str(l.split()[4])
		if l.split()[2]=="CDS":
			reg=re.search(r".*Parent=(.*)-R",l.split()[8])
			if reg.group(1) not in ghash:
				ghash[reg.group(1)+"%"+l.split()[6]]+=fastahash[l.split()[0]][int(l.split()[3])-1:int(l.split()[4])]
			if reg.group(1) in ghash:
				ghash[reg.group(1)+"_1%"+l.split()[6]]+=fastahash[l.split()[0]][int(l.split()[3])-1:int(l.split()[4])]
c=0	
clist=""
for k,v in ghash.items():
	
	if k.split("%")[1]=="+":
		genehash[k.split("%")[0]]=v
		for a in v:
			#print k, v[:12]
			if c!=3:
				clist+=a
				c+=1
			else:
				prothash[k.split("%")[0]]+=codon[clist]
				c=0
				clist=""
				clist+=a			
				c+=1
		prothash[k.split("%")[0]]+=codon[clist]
		clist=""
	if k.split("%")[1]=="_":
		genehash[k.split("%")[0]]=reversecomp(v)
		for a in reversecomp(v):
			#print k, v[:12]
			if c!=3:
				clist+=a
				c+=1
			else:
				prothash[k.split("%")[0]]+=codon[clist]
				c=0
				clist=""
				clist+=a			
				c+=1
		prothash[k.split("%")[0]]+=codon[clist]
		clist=""
out1=open(str(options.out)+"_prot.fasta","w")
c=0
clist=""
for k,v in prothash.items():
	
	out1.write(">"+infohash[k]+"\n")
	for a in v:
		if c!=60:
			clist+=a
			c+=1
		else:
			out1.write( clist+"\n")
			c=0
			clist=""
			clist+=a
			c+=1
	out1.write( clist+"\n")
	clist=""
	c=0
out2=open(str(options.out)+"_transcr.fasta","w")
c=0
clist=""
for k,v in genehash.items():
	
	out2.write( ">"+infohash[k]+"\n")
	for a in v:
		if c!=60:
			clist+=a
			c+=1
		else:
			out2.write( clist+"\n")
			c=0
			clist=""
			clist+=a
			c+=1
	out2.write( clist+"\n")
	clist=""
	c=0			
		