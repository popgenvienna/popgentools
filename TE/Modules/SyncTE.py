#!/usr/bin/env python
import sys
import collections


class TEPop:
	def __init__(self,presence,absence):
		self.__presence=presence
		self.__absence=absence
		
	@property
	def presence(self):
		return self.__presence
	
	@property
	def absence(self):
		return self.__absence
	
	@property
	def cov(self):
		return (self.__presence + self.__absence)

	@property
	def presfreq(self):
		if self.cov==0:
			return None
		return (float(self.presence)/float(self.cov))
		
	def __str__(self):
		return "{0}:{1}".format(str(self.presence),str(self.absence))
	


class SyncTE:
	
	def __init__(self,chrom, pos, support,family, order, knownid,comment,overlap,pops):
		self.__chrom=chrom
		self.__pos=pos
		self.__support=support
		self.__family=family
		self.__order=order
		self.__knownid=knownid
		self.__comment=comment
		self.__overlap=overlap
		self.__pops=pops

	
	@property
	def chrom(self):
		return self.__chrom
	
	@property
	def pos(self):
		return self.__pos
	
	@property
	def support(self):
		return self.__support
	
	@property	
	def family(self):
		return self.__family
	
	@property
	def order(self):
		return self.__order
	
	@property
	def knownid(self):
		return self.__knownid

	@property
	def comment(self):
		return self.__comment
	
	@property
	def overlap(self):
		return self.__overlap
	
	@property
	def pops(self):
		return self.__pops
	
	def __str__(self):
		toprint=[]
		toprint.append(self.chrom)
		toprint.append(str(self.pos))
		toprint.append(self.support)
		toprint.append(self.family)
		toprint.append(self.order)
		knownid=self.knownid
		if knownid is None:
			knownid="-"
		toprint.append(knownid)
		comment=self.comment
		if comment is None:
			comment="-"
		toprint.append(comment)
		toprint.append(str(self.overlap))
		for p in self.pops:
			toprint.append(str(p))
		return "\t".join(toprint)
		
	

		
	
# 4	149647	F	INE-1	TIR	FBti0062704	-	0	11:7	6:9	13:11
class SyncTEReader:
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")

	def __iter__(self):
		return self
	
	@classmethod
	def readfiltered(cls,file,ignoreoverlapping,mincoverage):
		teinsertions=SyncTEReader.readall(file)
		tes=[]
		for te in teinsertions:
			# discard overlapping ones
			if(te.overlap and ignoreoverlapping):
				continue
			
			# discard those below the minimum coverage
			mincovered=1
			for p in te.pops:
				if p.cov< mincoverage:
					mincovered=0
			if mincovered:
				tes.append(te)
		return tes
		
	@classmethod
	def readall(cls,file):
		sr=SyncTEReader(file)
		tes=[]
		for te in sr:
			tes.append(te)
		return tes

	
	def next(self):
		line=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":
				raise StopIteration
			line=line.rstrip('\n')
			if line != "" and line[0] != "#":
				break
		a=line.split()
		chrom=a.pop(0)
		pos=int(float(a.pop(0)))
		support=a.pop(0)
		family=a.pop(0)
		order=a.pop(0)
		knownid=a.pop(0)
		if knownid=="-":
			knownid=None
		comment=a.pop(0)	
		if comment=="-":
			comment=None
		overlap=int(a.pop(0))
		
		pops=[]
		for p in a:
			presence,absence=map(int,p.split(":"))
			po=TEPop(presence,absence)
			pops.append(po)
		
		return SyncTE(chrom, pos, support,family, order, knownid,comment,overlap,pops)