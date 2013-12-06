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
def loadhierarchy(file):
	"""
	KEPLER	FBgn0063570	Dbuz\Kepler	Foldback	Foldback	DNA
	1360	FBgn0005673	1360	SINE	non-LTR	RNA
	DNTOMRETA	FBgn0004357	Dana\Tom	LTR	LTR	RNA
	PPI251	FBgn0003055	P-element	SINE	non-LTR	RNA
	"""
	header=None
	hier={}
	leng=None
	for line in open(file):
		line=line.rstrip("\n")
		a=line.split("\t")
		alen=len(a)
		if leng is None:
			leng=alen
		if leng != alen:
			raise ValueError("All fields not have equal length in the categories")
		if(line.startswith("insert")):
			header=line
		else:
			a=line.split("\t")
			id=a.pop(0)
			if id in hier:
				raise ValueError("Invalid; " +id+" present multiple times")
			hier[id]=a
	return header,hier

def getTEsforChr(tean,chr):
	tesforchr=[]
	for te in tean:
		if te.chr == chr:
			tesforchr.append(te)
	return tesforchr

def get_uniquename(target,hier,prefix):
	counter=1
	name=target+"_"+prefix+"_m"+str(counter)
	while name in hier:
		counter+=1
		name=target+"_"+prefix+"_m"+str(counter)
	return name
		

def extractNovelSequences(noveltegtf,refseqs,hier,prefix):
	novelteseqs={}
	novelhier={}
	for chr in refseqs.keys():
		print "Processing chromosome "+chr+".."
		seq=refseqs[chr]
		tesforchr=getTEsforChr(noveltegtf,chr)
		for te in tesforchr:
			# process sequence
			start=te.start-1
			end=te.end
			novelseq=seq[start:end]			# get sequence of TE
			
			# update hierarchy
			target=te.target
			novelname=get_uniquename(target,hier,prefix)	# get a unique name
			novelteseqs[novelname]=novelseq		# add the sequence of the novel te to the database
			#print(novelname,target)
			hier[novelname]=hier[target] # add the novel hierarchy entry
			novelhier[novelname]=hier[target]
	return novelteseqs,novelhier


def printHierarch(header,hierarchy,outfile):
	ofw=open(outfile,"w")
	ofw.write(header+"\n")
	for k,v in hierarchy.items():
		a=[]
		a.append(k)
		a.extend(v)
		t="\t".join(a)
		ofw.write(t+"\n")
	ofw.close()

def printSequences(novelteseqs,outfasta):
	fw=FastaWriter(outfasta,100)
	idset=set([])
	for n,s in novelteseqs.items():
		if n in idset:
			raise ValueError("Already present"+n)
		fw.write(n,s)
	fw.close()
	
	

parser = OptionParser()
parser.add_option("--gtf-te",dest="gtfte",help="A gtf file containing the final TE annotation")
parser.add_option("--fasta-ref",dest="fastaref",help="A fasta file containing the reference sequence")
parser.add_option("--hierarchy",dest="hier",help="A hierarchy of TE sequences"),
parser.add_option("--out-hierarchy",dest="outhier",help="The output TE hierarchy"),
parser.add_option("--out-fasta",dest="outfasta",help="The output of the fasta sequences"),
parser.add_option("--prefix",dest="prefix",help="prefix of the novel TE sequences")
(options, args) = parser.parse_args()
prefix=options.prefix
if len(prefix)<1:
	raise ValueError("provide a prefix")
print("Loading refseqs..")
refseqs = FastaReader.readFastaHash(options.fastaref)
print("Loading gtf..")
noveltegtf= GTFTEReader.readall(options.gtfte)
print("Loading hierarchy..")
header,hierarchy= loadhierarchy(options.hier)

print("Extracting sequences and hierarchy..")
novelteseqs,novelhier=extractNovelSequences(noveltegtf,refseqs,hierarchy,prefix)

print("Printing hierarchy..")
printHierarch(header,novelhier,options.outhier)
print("Printing sequences..")
printSequences(novelteseqs,options.outfasta)
