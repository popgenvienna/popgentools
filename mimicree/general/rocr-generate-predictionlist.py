import sys
import random
from optparse import OptionParser, OptionGroup
import collections


class ResultReader:
	def __init__(self,file):
		self.__file=file
		self.__fh=open(file)
	
	"""
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
		majmin=SyncReader.parse_sync(a)
		return (chr,pos,majmin)
	"""
	def next(self):
		line=self.__fh.readline();
		if line=="":
			return None
		line=line.rstrip('\n');
		a=line.split()
		chr=a[0]
		pos=a[1]
		sig=a[-1]
		return (chr,pos,sig)

readar=[]
for i in range(1,len(sys.argv)):
	filename=sys.argv[i]
	readar.append(ResultReader(filename))


while(True):
	first=readar[0].next()
	if first is None:
		break
	activechr=first[0]
	activepos=first[1]
	topr=[first[2]]
	
	for i in range(1,len(readar)):
		tmp=readar[i].next()
		assert(activechr==tmp[0])
		assert(activepos==tmp[1])
		topr.append(tmp[2])
	toprstr="\t".join(topr)
	print(toprstr)
	
	
	
	
	