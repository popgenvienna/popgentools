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
		return (chr,int(pos),sig)
		
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
		return (chr,int(pos),sig)


class GuideWindowReader:
	def __init__(self,guidefile,windowsize):
		fh=None
		if(guidefile.endswith(".gz")):
			fh=gzip.open(guidefile,mode='rb')
		else:
			fh=open(guidefile)
		self.__gfh=fh
		self.__windowsize=windowsize
		self.__buffer=None
		self.__activechr=None
		self.__activepos=1
		
	def __iter__(self):
		return self
	
	def next(self):
		ac=self.__activechr
		startpos=self.__activepos
		endpos=startpos+self.__windowsize-1
		toret=[]
		while(1):
			line=self.__getnext()
			if line=="":
				raise StopIteration()
			line=line.rstrip("\n")
			a=line.split("\t")
			chr,pos=a[0],int(a[1])
			key=chr+":"+str(pos)
			if self.__activechr is None:
				self.__activechr=chr
				ac=chr
			if chr !=self.__activechr:
				self.__activechr=chr
				self.__activepos=1
				self.__bufferthis(line)
				break
			
			if pos <startpos:
				raise ValueError("invalid operation; file probably not sorted; {0} {1} {2}".format(pos,startpos,ac))
				
			if chr==self.__activechr and pos <=endpos:
				toret.append(key)
			elif(pos>endpos):
				self.__activepos=endpos+1
				self.__bufferthis(line)
				break
			else:
				raise ValueError("unhandled situation; should not occur")
			
		return (ac,startpos,toret)

	
	def __getnext(self):
		if self.__buffer is None:
			return self.__gfh.readline()
		else:
			toret=self.__buffer
			self.__buffer=None
			return toret
		
	def __bufferthis(self,sync):
		self.__buffer=sync



winsize=int(sys.argv[1])
chrguidefile=sys.argv[2]

readar=[]
for i in range(3,len(sys.argv)):
	filename=sys.argv[i]
	readar.append(ResultReader(filename))


for activechr,activepos,keys in GuideWindowReader(chrguidefile,winsize):
	if len(keys)==0:
		continue
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

	
	
	