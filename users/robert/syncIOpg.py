import sys
import random
import collections


class SyncWindowReader:
	def __init__(self,syncreader,windowsize):
		self.__sr=syncreader
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
			sync=self.__getnext()
			chr,pos=sync[0],sync[1]
			if self.__activechr is None:
				self.__activechr=chr
				ac=chr
			
			if chr !=self.__activechr:
				self.__activechr=chr
				self.__activepos=1
				self.__bufferthis(sync)
				break
			
			if pos <startpos:
				raise ValueError("invalid operation; file probably not sorted")
				
			if chr==self.__activechr and pos <=endpos:
				toret.append(sync)
			elif(pos>endpos):
				self.__activepos=endpos+1
				self.__bufferthis(sync)
				break
			else:
				raise ValueError("unhandled situation; should not occur")
			
		return (ac,startpos,endpos,toret)

	
	def __getnext(self):
		if self.__buffer is None:
			return self.__sr.next()
		else:
			toret=self.__buffer
			self.__buffer=None
			return toret
		
	def __buffer(self,sync):
		self.__buffer=sync

			
		

class SyncReaderMajMin:
	"""
	A light-weight sync-reader. Provides the population frequencies already as tuples of (majorallele, minorallele)
	
	returns chromosome, position, (major-allele, minor-allele), ((pop1-maj,pop1-min),(pop2-maj,pop2-min),...)
	
	For example
	2L	15	A	45:0:108:0:0:0	47:0:90:0:0:0
	produces
	[2L,15,C,A),((108,45),(90,47))]
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
		pos=int(a.pop(0))
		refc=a.pop(0)
		mami,majmin=SyncReaderMajMin.__parse_sync(a)
		return (chr,pos,mami,majmin)

		
	@classmethod
	def __parse_sync(cls,entries):
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
		return [(major,minor),toret]

class SyncReaderRevAllele:
	"""
	A light-weight sync-reader. Provides the count of the reference allele and the coverage
	
	returns chromosome, position, refChar ((rcount_p1,cov_p1), (rcount_p2,cov_p2),...)
	
	For example
	2L	15	A	45:0:108:0:0:0	47:0:90:0:0:0
	produces
	[2L, 15, A, ((5,10),(10,12))]
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
		pos=int(a.pop(0))
		refc=a.pop(0)
		refc=refc.upper() #convert to upper case
		reff=self.__parse_sync_af(a,refc)
		return (chr,pos,refc,reff)

		
	
	def __parse_sync_af(self,entries,refc):
		parsed=[]
		for e in entries:
			a=map(float,e.split(":")) # float conversions
			np={'A':a[0],'T':a[1],'C':a[2],'G':a[3]}
			cov=a[0]+a[1]+a[2]+a[3]

			if(refc not in np):
				parsed.append((0,cov))
			else:
				parsed.append((np[refc],cov))
		return parsed


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

class SyncWriterATCG:
	"""
	A light-weight sync writer; requires the counts in ATCGNDel
	requires chromosome, position, refChar ([A,T,C,G,N,del], [A,T,C,G,N,del],...) # first ATCGNdel for first population, second for second ...
	"""
	def __init__(self,file):
		self.__filename=file
		sefl.__filehandle=open(file,"w")
	
	def write(self,chr,pos,refChar,pops):
		topr=SyncWriterATCG.format(chr,pos,refChar,pops)
		self.__filehandle.write(topr+"\n")
	
	@classmethod
	def format(cls,chr,pos,refChar,pops):
		tojoin=[str(chr),str(pos),str(refChar)]
		formpops=[]
		for p in pops:
			psi=[str(int(i)) for i in p]
			fp=":".join(psi)
			formpops.append(fp)
		tojoin.extend(formpops)
		topr="\t".join(tojoin)
		return topr
		
		
	def close():
		self.__filehandle.close()
