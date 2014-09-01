import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import gzip


class ResultReader:
	def __init__(self,file):
		self.__file=file
		fh=None
		if(file.endswith(".gz")):
			fh=gzip.open(file,mode='rb')
		else:
			fh=open(file)
		self.__fh=fh
		self.__buffer=""
		self.__lastline=""
		
	def __iter__(self):
		return self
	
	def __getnext(self):
		if(self.__buffer!=""):
			toret=self.__buffer
			self.__buffer=""
			return toret
		else:
			return self.__fh.readline()
			
			
	def buffer(self):
		self.__buffer=self.__lastline
		
	def read(self):
		line=self.__getnext();
		if line=="":
			return ("","","")
		self.__lastline=line
		line=line.rstrip('\n');
		a=line.split()
		chr=a[0]
		pos=a[1]
		sig=a[-1]
		return (chr,pos,sig)
		
	def next(self):
		line=self.__getnext();
		if line=="":
			raise StopIteration
		self.__lastline=line
		line=line.rstrip('\n');
		a=line.split()
		chr=a[0]
		pos=a[1]
		sig=a[-1]
		return (chr,pos,sig)

guidereader=ResultReader(sys.argv[1])
readar=[]
for i in range(2,len(sys.argv)):
	filename=sys.argv[i]
	readar.append(ResultReader(filename))


for activechr,activepos,sig in guidereader:
	topr=[]
	for i in range(0,len(readar)):
		chr,pos,sig=readar[i].read()
		if(chr==activechr and pos==activepos):
			topr.append(sig)
		else:
			topr.append("0.0")
			readar[i].buffer()
	toprstr="\t".join(topr)
	print(toprstr)
	
	
	
	
	