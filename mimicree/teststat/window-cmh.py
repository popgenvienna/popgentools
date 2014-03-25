#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections


class GuideWindowCMHReader:
	def __init__(self, guidefile, windowsize):
		self.__gfh=open(guidefile)
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
			key=float(a.pop())
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

			



def parse_line(syncs,baselist,derivedlist):
	# (C,A),(108,45),(90,47))
	dera=[]
	basea=[]
	for chr,pos,mami,parsed in syncs:
		derc=[0,0]
		basec=[0,0]
		for i in derivedlist:
			key=i-1
			derc[0]+=parsed[key][0]
			derc[1]+=parsed[key][1]
		for i in baselist:
			key=i-1
			basec[0]+=parsed[key][0]
			basec[1]+=parsed[key][1]
		dera.append(derc)
		basea.append(basec)
	return (basea,dera)




		

parser = OptionParser()
parser.add_option("--sync",dest="sync",help="A file containing the cmh results")
parser.add_option("--compare",dest="comp",help="A comparision string as with the cmh-test")
(options, args) = parser.parse_args()

for chr,start,cmhs in GuideWindowCMHReader(options.sync,1000):
	if len(cmhs)==0:
		continue
	sum=0.0
	for c in cmhs:
		sum+=c
	avcmh=sum/float(len(cmhs))
	
	topr=[chr,str(start),str(avcmh)]
	print "\t".join(topr)
