import sys
import random
from optparse import OptionParser, OptionGroup
import collections


class SyncReaderATCG:
	"""
	A light-weight sync reader; Provides the counts of ATCGNdel
	
	returns chromosome, position, refChar ([A,T,C,G,N,del], [A,T,C,G,N,del],...) # first ATCGNdel for first population, second for second ...
	"""
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")
	

	def __iter__(self):
		return self
	
	def next(self):
		line=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":
				raise StopIteration
			line=line.rstrip('\n')
			if line != "":
				break
		
		a=line.split()
		chr=a.pop(0)
		pos=a.pop(0)
		refc=a.pop(0)
		counts=SyncReaderATCG.parse_sync_cov(a)
		return (chr,pos,refc,counts)
		
	@classmethod
	def parse_sync_cov(cls,entries):
		conv=[]
		for e in entries:
			a=map(float,e.split(":"))
			conv.append(a)
		return conv


parser = OptionParser()
parser.add_option("--input",dest="input",help="the input file as sync")
(options, args) = parser.parse_args()

covs=collections.defaultdict(lambda:0)

for t1,t2,t3,entries in SyncReaderATCG(options.input):
	for e in entries:
		cov=e[0]+e[1]+e[2]+e[3]
		cov=int(cov)
		covs[cov]+=1

for k,v in covs.items():
	print "{0}\t{1}".format(k,v)

