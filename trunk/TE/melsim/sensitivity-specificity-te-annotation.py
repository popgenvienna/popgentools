from gtfIO import GTFReader,GTFEntry;
import sys
import random
from optparse import OptionParser, OptionGroup
import collections

def convert_chrhash(entries):
	toret=collections.defaultdict(lambda:[])
	for e in entries:
		toret[e.chr].append(e)
	return toret

def create_basehash(entries):
	bh={}
	for e in entries:
		for i in range(e.start,e.end+1):
			bh[i]=1
	return bh

def senspec_base(cref,ctotest):
	bhref=create_basehash(cref)
	bhtotest=create_basehash(ctotest)
	countcommon=0
	for k in bhref.keys():
		if k in bhtotest:
			countcommon+=1
	countbhref=len(bhref.keys())
	countbhtotest=len(bhtotest.keys())
	return countcommon,countbhref,countbhtotest

def countBinA(a,b):
	sref=sorted(a,key=lambda te: te.start)  # a is the reference
	# b does not need to be sorted; irrelevant
	countoverlap=0
	countnotoverlap=0
	for tt in b:
		# testing every te in b (totest)
		overlap=False
		for ref in sref:
			# STOP Condition: start position of the reference is already larger than the end position of the tested TE
			if ref.start> tt.end:
				break
			if((ref.start <= tt.start and ref.end >= tt.start) or (ref.start <= tt.end and ref.end >= tt.end)):
				overlap=True
				break
		if(overlap):
			countoverlap+=1
		else:
			countnotoverlap+=1
	return (countoverlap,countnotoverlap)

def identifyOverlap(a):
	sref=sorted(a,key=lambda te: te.start)
	countoverlap=0
	for i in range(0,len(sref)):
		for k in range(i+1,len(sref)):
			tt=sref[i]
			ref=sref[k]
			if ref.start> tt.end:
				break
			if((ref.start <= tt.start and ref.end >= tt.start) or (ref.start <= tt.end and ref.end >= tt.end) or (ref.start > tt.start and ref.end < tt.end)):
				countoverlap+=1
				break
	return countoverlap

def senspec_entries(cref,ctotest):
	# >>> sorted(student_tuples, key=lambda student: student[2])   # sort by age

	coverlap,cnottotest=countBinA(cref,ctotest)
	coverlap2,cnotref=countBinA(ctotest,cref)
	
	#assert(coverlap<=coverlap2) # shity reference annotation may contain overlapping TE annotation...
	return(coverlap,coverlap2,cnotref,cnottotest)

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
			target=target.lstrip('"Motif:')
			target=target.rstrip('"')
			e.target=target
			return e
		raise StopIteration

def filterregion(tofilter,region):
	filtered=[]
	for te in tofilter:
		if te.chr not in region:
			continue
		r=region[te.chr]
		if(te.start >= r[0] and te.end <= r[1]):
			filtered.append(te)
	return filtered


parser = OptionParser()
parser.add_option("--reference",dest="reference",help="A gtf file containing the reference annotation")
parser.add_option("--totest",dest="totest",help="A gtf file containing the annotation that should be tested")
(options, args) = parser.parse_args()
regionfilter={
	"X":  [0,22422827],
	"2L": [0,22420241], 
	"2R": [387345,21146708],
	"3L": [0,23825333],
	"3R": [0,27905053],
	"4" : [0,1350078]}
# chromosomes with incomplete reference annotation: 2L, 2R, 3L
# fine: 3R, X, 4


rawref=filterregion(GTFReader.readall(options.reference),regionfilter)
rawtotest=filterregion(GTFReader.readall(options.totest),regionfilter)
chref=convert_chrhash(rawref)
chtotest=convert_chrhash(rawtotest)

# ONLY TEST the chromosomes in the given list
tcboverlap,tcbref,tcbtotest=(0,0,0) # total count bases =tcb
tceoverlap,tceoverlap2,tceref,tcetotest=(0,0,0,0) # total count element =tcb
countRef=0
countTotest=0
countBRef=0
countBTotest=0
for chr in regionfilter.keys():
	cref=chref[chr]
	ctotest=chtotest[chr]
	#intoverlapref=identifyOverlap(cref)
	#intoverlaptotest=identifyOverlap(ctotest)
	#print(chr,intoverlapref,intoverlaptotest)
	countRef+=len(cref)
	countTotest+=len(ctotest)
	eoverlap,eoverlap2,eref,etotest=senspec_entries(cref,ctotest)
	tceoverlap+=eoverlap
	tceoverlap2+=eoverlap2
	tceref+=eref
	tcetotest+=etotest
	boverlap,breftot,btotesttot=senspec_base(cref,ctotest)
	countBRef+=breftot
	countBTotest+=btotesttot
	bref=breftot-boverlap
	btotest=btotesttot-boverlap
	tcboverlap+=boverlap
	tcbref+=bref
	tcbtotest+=btotest
	

print "{0}\t{1}\t{2}\t{3}".format(countRef,countTotest,countBRef,countBTotest)
print "{0}\t{1}\t{2}".format(tcboverlap,tcbref,tcbtotest)
print "{0}\t{1}\t{2}".format(tceoverlap2,tceref,tcetotest)


	


	

	
	



