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

        

def parse_line(entries,baselist,derivedlist):
        parsed=entries
        
        derivedentries=[]
        baseentries=[]
        for i in derivedlist:
                key=i-1
                derivedentries.append(parsed[key])
        
        for i in baselist:
                key=i-1
                baseentries.append(parsed[key])
        return (baseentries,derivedentries)
        
def get_associastat(base,derived):
        assert(len(base)==len(derived))
        
        associastat=0
        for i,b in enumerate(base):
                d=derived[i]
                bf=b[0]/(b[0]+b[1])
                df=d[0]/(d[0]+d[1])
                diff=df-bf
                associastat+=diff
        associastat=abs(associastat)
        return associastat

def parse_comparestring(comp):
	a=comp.split(",")
	
	baselist=[]
	derivedlist=[]
	for e in a:
		b,d=map(int,e.split("-"))
		baselist.append(b)
		derivedlist.append(d)
	return (baselist,derivedlist)
		

parser = OptionParser()
parser.add_option("--sync",dest="sync",help="A file containing the cmh results")
parser.add_option("--compare",dest="comp",help="A comparision string as with the cmh-test")
(options, args) = parser.parse_args()

baselist,derivedlist=parse_comparestring(options.comp)

for chr,pos,s in SyncReader(options.sync):
        baseentries,derivedentries=parse_line(s,baselist,derivedlist)
        assstat=get_associastat(baseentries,derivedentries)
        
        topr=[chr,pos,str(assstat)]
        print "\t".join(topr)
        