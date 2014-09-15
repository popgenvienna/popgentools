import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import gzip




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

			



def get_selected_hash(file):
	s=set([])
	for line in open(file):
		"""
		2L	6009099	A	0.1	0.5
		3R	12520973	A	0.1	0.5
		"""
		line=line.rstrip()
		a=line.split("\t")
		chr,pos=(a[0],a[1])
		key=chr+":"+pos
		s.add(key)
	return s


chrguidefile=sys.argv[1]

selectedar=[] # content will be sets with the key: "chr:pos" as string
for i in range(2,len(sys.argv)):
	filename=sys.argv[i]
	selectedar.append(get_selected_hash(filename))


for chr,pos,keys in GuideWindowReader(chrguidefile,1000):
	if len(keys)==0:
		continue
	
	topr=[]
	for s in selectedar:
		topr.append(0);
	
	for key in keys:
		for i,s in enumerate(selectedar):
			if key in s:
				topr[i]=1
	toprstr="\t".join(map(str,topr))
	print toprstr

	
	
	
	