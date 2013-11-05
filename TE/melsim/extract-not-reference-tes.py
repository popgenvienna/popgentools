from gtfIO import GTFReader,GTFEntry,GTFWriter;
import sys
import random
from optparse import OptionParser, OptionGroup
import collections

def convert_chrhash(entries):
	toret=collections.defaultdict(lambda:[])
	for e in entries:
		toret[e.chr].append(e)
	return toret


def getBinA(a,b):
	sref=sorted(a,key=lambda te: te.start)  # a is the reference
	# b does not need to be sorted; irrelevant
	
	novel=[]
	for tt in b:
		# testing every te in b (totest)
		overlap=False
		for ref in sref:
			# STOP Condition: start position of the reference is already larger than the end position of the tested TE
			if ref.start> tt.end:
				break
			if((ref.start <= tt.start and ref.end >= tt.start) or (ref.start <= tt.end and ref.end >= tt.end) or (ref.start > tt.start and ref.end < tt.end)):
				overlap=True
				break
		if( not overlap):
			novel.append(tt)

	return (novel)


def novelentries(cref,ctotest):
	# >>> sorted(student_tuples, key=lambda student: student[2])   # sort by age
	novel=getBinA(cref,ctotest)
	return(novel)

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


parser = OptionParser()
parser.add_option("--reference",dest="reference",help="A gtf file containing the reference annotation")
parser.add_option("--totest",dest="totest",help="A gtf file containing the annotation that should be tested")
parser.add_option("--output",dest="output",help="The output file containing novel TE insertions not in the reference genome")
(options, args) = parser.parse_args()

# chromosomes with incomplete reference annotation: 2L, 2R, 3L
# fine: 3R, X, 4


rawref=GTFReader.readall(options.reference)
rawtotest=GTFReader.readall(options.totest)
chref=convert_chrhash(rawref)
chtotest=convert_chrhash(rawtotest)

ofh=GTFWriter(options.output)

for chr in ["X","2L","2R","3L","3R","4"]:
	cref=chref[chr]
	ctotest=chtotest[chr]
	ne=novelentries(cref,ctotest)
	for n in ne:
		ofh.write(n)
ofh.close()


	


	

	
	



