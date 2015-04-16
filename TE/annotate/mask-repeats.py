from gtfIO import GTFReader,GTFEntry;
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from fastaIO import FastaReader,FastaWriter

class GTFTEReader:
	def __init__(self,file):
		self.__gtfreader=GTFReader(file)
	
	def __iter__(self):
		return self
	
	def next(self):
		for e in self.__gtfreader:
			# 	2L	RepeatMasker	similarity	1724	1849	 6.4	-	.	Target "Motif:DMTRDNA" 1311 1435
			a=e.comment.split(" ") 
			target=a[1]
			score=float(e.end-e.start)
			e.sore=score
			target=target.replace('"Motif:',"")
			target=target.rstrip('"')
			e.target=target
			return e
		raise StopIteration
	
	@classmethod
	def readall(cls,file):
		toret=[]
		for e in GTFTEReader(file):
			toret.append(e)
		return toret

def maskSeq(seq,start,end):
	first=seq[0:start]
	last=seq[end:]
	middle="N"*(end-start)
	return first+middle+last

def getTEsforChr(tean,chr):
	tesforchr=[]
	for te in tean:
		if te.chr == chr:
			tesforchr.append(te)
	return tesforchr

def get_uniquename(target,hier):
	counter=1
	name=target+"_m"+str(counter)
	while name in hier:
		counter+=1
		name = target+"_m"+str(counter)
	return name
		

def maskTEsinSeq(tean,refseq):
	novelrefseq={}
	for chr in refseq.keys():
		print "Processing chromosome "+chr+".."
		seq=refseq[chr]
		tesforchr=getTEsforChr(tean,chr)
		for te in tesforchr:
			# process sequence
			start=te.start-1
			end=te.end
			seq=maskSeq(seq,start,end)		# mask the TE with Ns
		novelrefseq[chr]=seq
	return novelrefseq



def printSequences(seq,outfasta):
	fw=FastaWriter(outfasta,60)
	for n,s in seq.items():
		fw.write(n,s)
	fw.close()
	
	

parser = OptionParser()
parser.add_option("--gtf",dest="gtfte",help="A gtf file containing the TE annotation")
parser.add_option("--input",dest="fastaref",help="A fasta file containing the reference sequence")
parser.add_option("--output",dest="outfasta",help="The output of the fasta sequences"),


(options, args) = parser.parse_args()
print("Loading refseqs..")
refseqs = FastaReader.readFastaHash(options.fastaref)
print("Loading gtf..")
noveltegtf= GTFTEReader.readall(options.gtfte)
print("Masking reference sequence..")
novelrefseq=maskTEsinSeq(noveltegtf,refseqs)
print("Printing masked reference sequence..")
printSequences(novelrefseq,options.outfasta)
