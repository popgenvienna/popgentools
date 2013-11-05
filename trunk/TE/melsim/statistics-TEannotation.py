from gtfIO import GTFReader,GTFEntry;
import sys
import random
from optparse import OptionParser, OptionGroup
import collections

parser = OptionParser()
parser.add_option("--gtf",dest="gtf",help="A gtf file containing the RepeatMasked gtf annotation")
parser.add_option("--id",dest="id",help="Sample ID")
(options, args) = parser.parse_args()



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

id=options.id
count=0
lengthsum=0
lengdistri=collections.defaultdict(lambda:0)
famhash=collections.defaultdict(lambda:0)


for e in GTFTEReader(options.gtf):
	count+=1
	leng=e.end-e.start+1
	lengthsum+=leng
	famhash[e.target]+=1
	lengdistri[leng]+=1

print "Count\t{0}".format(count)
print "Total TE \t{0}".format(lengthsum)
print "Average length {0}".format(float(lengthsum)/float(count))
famitems=sorted(famhash.items(),key=lambda k:-k[1])
for fi in famitems:
	print "{0} {1}".format(fi[0],fi[1])

for ld in sorted(lengdistri.keys()):
	print "LDi\t{0}\t{1}\t{2}".format(id,ld,lengdistri[ld])

	




