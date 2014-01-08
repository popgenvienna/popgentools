import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import re
from gtfIO import GTFEntry, GTFWriter


class NucmerEntry:
	def __init__(self,rname,qname,rstart,rend,qstart,qend):
		self.rname=rname
		self.qname=qname
		if(rend < rstart):
			(rstart,rend)=(rend,rstart)
		if(qend < qstart):
			(qstart,qend)=(qend,qstart)
		self.rstart=rstart
		self.rend=rend
		self.qstart=qstart
		self.qend=qend
		self.dist=[]
		self.count=1

def filter_chr(tofilter,chrnames):
	filtered=[]
	for f in tofilter:
		if f.rname in chrnames and f.qname in chrnames:
			filtered.append(f)
	return filtered

def get_chrhash(tohash):
	chrh=collections.defaultdict(lambda:[])
	for t in tohash:
		chrh[t.rname].append(t)
	return chrh
			
			
def get_distance(nm1,nm2):
	"""
	rname,qname,rstart,rend,qstart,qend
	"""
	distunqual=1000000000
	if(nm1.rname != nm2.rname):
		return distunqual, distunqual
	if(nm1.qname != nm2.qname):
		return distunqual, distunqual
	distr = _getdist(nm1.rstart, nm1.rend, nm2.rstart, nm2.rend)
	distq = _getdist(nm1.qstart, nm1.qend, nm2.qstart, nm2.qend)
	if(distq>distr):
		return distq,distr
	else:
		return distr,distr
	
def _getdist(s1,e1,s2,e2):
	if(s2>e1):
		return s2-e1
	elif(e2<s1):
		return s1-e2
	else:
		return 0

def mergeNucmer(nm1,nm2,dist):
	assert nm1.rname==nm2.rname
	assert nm1.qname==nm2.qname
	
	(nue_rstart,nue_rend,nue_qstart,nue_qend) = (nm1.rstart,nm1.rend,nm1.qstart,nm1.qend)
	if(nm2.rstart<nue_rstart):
		nue_rstart=nm2.rstart
	if(nm2.rend>nue_rend):
		nue_rend=nm2.rend
	if(nm2.qstart<nue_qstart):
		nue_qstart=nm2.qstart
	if(nm2.qend>nue_qend):
		nue_qend=nm2.qend
	nuevo=NucmerEntry(nm1.rname,nm1.qname,nue_rstart,nue_rend,nue_qstart,nue_qend)
	nuevo.dist.extend(nm1.dist)
	nuevo.dist.extend(nm2.dist)
	nuevo.dist.append(dist)
	nuevo.count=nm1.count+nm2.count
	return nuevo

def get_refgtflist(nucmerentries):
	#chr, source, feature, start, end, score, strand, frame, comment (comment is unparsed)
	gtfes=[]
	
	for nm in nucmerentries:
		comment="{0}:{1}-{2} count={3}".format(nm.qname,nm.qstart,nm.qend,nm.count)
		gtfes.append(GTFEntry(nm.rname,"nucmer","orthologous",nm.rstart,nm.rend, ".", ".", ".", comment))
	return gtfes

def get_querygtflist(nucmerentries):
	#chr, source, feature, start, end, score, strand, frame, comment (comment is unparsed)
	gtfes=[]
	for nm in nucmerentries:
		comment="{0}:{1}-{2} count={3}".format(nm.rname,nm.rstart,nm.rend,nm.count)
		gtfes.append(GTFEntry(nm.qname,"nucmer","orthologous",nm.qstart,nm.qend, ".", ".", ".", comment))
	return gtfes

def print_distances(file,toprint):
	ofw=open(file,"w")
	
	for tp in toprint:
		for d in tp.dist:	
			ofw.write("{0}\n".format(d))
	ofw.close()

def read_nucmer(input):
	"""
	[S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
	=====================================================================================
	181689   181893  |       55      251  |      205      197  |    88.78  | U	contig_1116
	18556    18794  |        1      253  |      239      253  |    90.62  | Uextra	contig_3214
	17790557 17790688  |     1678     1809  |      132      132  |    96.21  | Uextra	contig_3214
	
	From the manual
	Output is to stdout and is slightly different depending on the type of alignment, i.e. nucleotide or amino acid.
	Some of the described columns, such as percent similarity, will not appear for nucleotide comparisons.
	When run without the -H or -B options, show-coords prints a header tag for each column;
	the descriptions of each tag follows.
	[S1] start of the alignment region in the reference sequence
	[E1] end of the alignment region in the reference sequence
	[S2] start of the alignment region in the query sequence
	[E2] end of the alignment region in the query sequence
	[LEN 1] length of the alignment region in the reference sequence [LEN 2] length of the alignment region in the query sequence[% IDY] percent identity of the alignment [% SIM] percent similarity of the alignment (as determined by the BLOSUM scoring matrix)[% STP] percent of stop codons in the alignment [LEN R] length of the reference sequence [LEN Q] length of the query sequence[COV R] percent alignment coverage in the reference sequence [COV Q] percent alignment coverage in the query sequence[FRM] reading frame for the reference and query sequence alignments respectively
	[TAGS] the reference and query FastA IDs respectively. All output coordinates and lengths are relative to the forward strand of the reference DNA sequence.
	"""
	readflag=0
	toret=[]
	for line in open(input):
		line=line.rstrip("\n")
		if line.startswith("====================================================================================="):
			readflag=1
		elif(readflag):
			a=line.split('|')
			tmp=a[0].strip()
			refcoords=re.split(r"\s+",tmp)
			tmp=a[1].strip()
			querycoords=re.split(r"\s+",tmp)
			tmp=a[4].strip()
			names=re.split(r"\s+",tmp)
			nue=NucmerEntry(names[0],names[1],int(refcoords[0]),int(refcoords[1]),int(querycoords[0]),int(querycoords[1]))
			toret.append(nue)
	return toret

parser = OptionParser()
parser.add_option("--nm",dest="nucmer",help="A nucmer show coordindates file")
parser.add_option("--output-ref",dest="outputref",help="the output file")
parser.add_option("--max-dist",dest="maxdist",help="the  maximum distance between two fragments")
parser.add_option("--output-query",dest="outputquery",help="the output query")
parser.add_option("--output-distances",dest="outputdist",help="the output of the distances")
(options, args) = parser.parse_args()
euchromosomes=set(["X","2L","2R","3L","3R","4"])

maxdist=int(options.maxdist)
rawnucmer=read_nucmer(options.nucmer)
print "Read {0} alignments".format(len(rawnucmer))
filternucmer=filter_chr(rawnucmer,euchromosomes)
print "Retained {0} alignments on euchromosomes".format(len(filternucmer))

finallist=[]
chrhash=get_chrhash(filternucmer)
comparecount=1
for chr, chrlist in chrhash.items():
	
	print "Processing {0} having {1} fragments".format(chr,len(chrlist))
	chrlist=sorted(chrlist,key=lambda k:k.rstart)
	chrentries=[]
	while(len(chrlist)>0):
		active=chrlist.pop(0)
		
		while(True):
			tested=[]
			activeUnmodified=True
			while(len(chrlist)>0):
				totest=chrlist.pop(0)
				#if(comparecount%100000==0):
				#	print "Finished {0} comparision; len chrlist {1}".format(comparecount,len(chrlist))
				comparecount+=1
				testdist,narrowdist=get_distance(active,totest)
				if(testdist < maxdist):
					active = mergeNucmer(active,totest,testdist)
					activeUnmodified=False
				else:
					tested.append(totest)
			chrlist=tested
			if activeUnmodified:
				break
		chrentries.append(active)
	print "Finished {0} with {1} fragments".format(chr,len(chrentries))
	finallist.extend(chrentries)
refgtf=get_refgtflist(finallist)
querygtf=get_querygtflist(finallist)
GTFWriter.write_all(options.outputref,refgtf) # (cls,file,gtfentries):
GTFWriter.write_all(options.outputquery,querygtf)

if(options.outputdist):
	print_distances(options.outputdist,finallist)

