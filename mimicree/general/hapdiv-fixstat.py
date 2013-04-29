#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections

class SyncReader:
	"""
	Light weight sync reader
	2L	15	A	45:0:108:0:0:0	47:0:90:0:0:0	52:0:103:0:0:0	58:0:91:0:0:0	57:0:107:0:0:0	49:0:84:0:0:0
	2L	47	A	155:0:0:0:0:0	141:0:0:0:0:0	152:0:0:0:0:0	146:0:0:0:0:0	161:0:1:0:0:0	172:0:0:0:0:0
	2L	57	A	148:0:4:0:0:0	128:0:4:0:0:0	172:0:4:0:0:0	164:0:1:0:0:0	148:0:4:0:0:0	141:0:1:0:0:0
	2L	89	A	156:0:3:0:0:0	161:0:6:0:0:0	146:0:5:0:0:0	158:0:4:0:0:0	163:0:6:0:0:0	151:0:0:0:0:0
	2L	224	A	154:0:0:0:0:0	149:0:0:0:0:0	135:0:0:0:0:0	167:0:0:0:0:0	164:0:0:0:0:0	162:0:0:0:0:0
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
		majmin=SyncReader.parse_sync(a)
		return (chr,pos,majmin)

		
	@classmethod
	def parse_sync(cls,entries):
		parsed=[]
		for e in entries:
			a=map(float,e.split(":"))
			np={'A':a[0],'T':a[1],'C':a[2],'G':a[3]}
			parsed.append(np)
		ac,tc,cc,gc = (0,0,0,0)
		
		for p in parsed:
			ac += p['A']
			tc += p['T']
			cc += p['C']
			gc += p['G']
			
		tmpar=[ (ac,'A'),
			(tc,'T'),
			(cc,'C'),
			(gc,'G') ]
		
		tmpar=sorted(tmpar, key=lambda cs: -cs[0])
		major=tmpar[0][1]
		minor=tmpar[1][1]
		toret=[]
		for p in parsed:
			novel=(p[major],p[minor])
			toret.append(novel)
		return toret


        



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

parser = OptionParser()
selected=sys.argv[1]

selhash=get_selected_hash(selected)


filename=sys.argv[2]
selected=0
fixed=0
for line in open(filename):
        line=line.rstrip()
        a=line.split("\t")
        chr,pos=(a[0],a[1])
        key=chr+":"+pos
        if key in selhash:
                selected+=1
                parsed=SyncReader.parse_sync([a[3]])[0]
                if(parsed[0]==0 or parsed[1]==0):
                        fixed+=1
topr=[]
topr.append(filename)
topr.append(str(selected))
topr.append(str(fixed))
toprstr="\t".join(topr)
print(toprstr)
