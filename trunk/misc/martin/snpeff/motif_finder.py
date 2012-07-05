from Bio import SeqIO
import sys
import collections



def fasta_parse(x):
	h=collections.defaultdict(lambda:"")	
	genome=SeqIO.parse(open(x), "fasta")
	for record in genome:
		h[record.id]=record.seq

	return h

f=fasta_parse(sys.argv[1])
features=["INTRON","UPSTREAM","DOWNSTREAM","UTR_5_PRIME","UTR_3_PRIME","INTERGENIC"]
for l in open(sys.argv[2],"r"):
	a=l.split("\t")
	if a[-4] in features:
		print ">"+a[0]+str(a[2])+"_"+a[-4]
		print f[a[0]][int(a[2])-11:int(a[2])+9]