#!/usr/bin/env python


class MimicreeHap:
	"""
	Represents a line in a haplotype file, which corresponds to the genotypes of a whole population
	X       9121094 T       T/C     TC CT CC CT CT CC TT CT TT TC CC CC CC CT TT
	"""
	def __init__(self,chr,pos,refc,major,minor,haplotypes):	
		self.chr=chr
		self.pos=pos
		self.refc=refc
		self.major=major
		self.minor=minor
		self.__haplotypes=haplotypes
		self.__genotypes=None

	
	@property
	def haplotypes(self):
		return self.__haplotypes
	
	@property
	def genotypes(self):
		if(self.__genotypes is None):
			haps=self.__haplotypes
			genotypes=[haps[i]+haps[i+1] for i in range(0,len(haps),2)] # nice despite slight hangover ;)
			self.__genotypes=genotypes
		return self.__genotypes
	
	@property
	def haplotypeCount(self):
		return len(self.haplotypes)
		
	@property
	def genotypeCount(self):
		return len(self.genotypes)
		
	@property
	def countMajMin(self):
		maj,min=(0,0)
		for i in  self.haplotypes:
			if(i ==self.major):
				maj+=1
			elif(i==self.minor):
				min+=1
			else:
				raise ValueError("Invalid allele")
		return (maj,min)
	
		
	def countAllele(self,allele):
		a=0
		for i in self.haplotypes:
			if(i==allele):
				a+=1
		return a
	
	def freqAllele(self,allele):
		ca=float(self.countAllele(allele))
		sum=float(self.haplotypeCount)
		return ca/sum

	@property
	def freqMajMin(self):
		cmaj,cmin=self.countMajMin
		sum=float(cmaj+cmin)
		fmaj=float(cmaj)/sum
		fmin=float(cmin)/sum
		return (fmaj,fmin)
	
	@property
	def isPolymorphic(self):
		alleleset=set(self.haplotypes)
		assert(len(alleleset)>0)
		if(len(alleleset)>1):
			return True
		else:
			return False

class MimicreeHapWindowReader:
	"""
	Reads a MimicrEE Haplotype file in chunks of the given window size
	returns a tuple consisting of chr, pos, window
	chr.. chromosome
	pos.. position
	window.. all entries for the given window as class MimicreeHap
	"""
	
	def __init__(self,file,windowsize):
		self.__windowsize=windowsize
		self.__mhr=MimicreeHapReader(file)
		self.__buffer=[]
		self.__activeChr=None
		self.__lastPos=windowsize
	
	def __iter__(self):
		return self
		
	def __getNext(self):
		if(len(self.__buffer)>0):
			return self.__buffer.pop(0)
		n=None
		try:
			n=self.__mhr.next()
		except StopIteration:
			pass
		return n
	
		
	def __pushBuffer(self,tobuffer):
		self.__buffer.append(tobuffer)
	
	
	def next(self):
		toret=[]
		toretchr=self.__activeChr
		toretpos=self.__lastPos - (float(self.__windowsize)/2.0)
		while(1):
			n=self.__getNext()
			if(n is None):
				break
			if(self.__activeChr is None):
				self.__activeChr = n.chr
				toretchr=n.chr
			
			if(n.chr != self.__activeChr):
				self.__pushBuffer(n)
				self.__lastPos=self.__windowsize
				self.__activeChr=n.chr
				return toretchr,toretpos,toret
			if(n.pos>self.__lastPos):
				self.__pushBuffer(n)
				self.__lastPos+=self.__windowsize
				return toretchr,toretpos,toret
			
			toret.append(n)
	
		if(len(toret)>0):
			return toretchr, toretpos, toret
		raise StopIteration


class MimicreeHapReader:
	
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
		
		line=line.rstrip()
		# the external parse method may be quite useful in many circumstances, so do not touch
		toret = MimicreeHapIO.parseLine(line)
		return toret
	
	

class MimicreeHapIO:
	
	@classmethod
	def genotypeCount(cls,inputFile):
		"""
		Return the number of genotypes present in a input file
		"""
		fh=open(inputFile)
		line=fh.next()
		fh.close()
		pa=MimicreeHapIO.parseLine(line)
		return pa.genotypeCount
	
	@classmethod
	def haplotypeCount(cls,inputFile):
		fh=open(inputFile)
		line=fh.next()
		fh.close()
		pa = MimicreeHapIO.parseLine(line)
		return pa.haplotypeCount
	

	@classmethod
	def parseLine(cls,line):
		a=line.split("\t");
		(major,minor)=a[3].split("/")
		genotypes=a[4].split(" ")
		haplotypes=[]
		for g in genotypes:
			haplotypes.append(g[0])
			haplotypes.append(g[1])
		#(lengen,lenhap)=len(genotypes),len(haplotypes)
		#test=genotypes[-1]
		#test2=haplotypes[-1]
		#test3=haplotypes[-2]
		return MimicreeHap(a[0],int(a[1]),a[2],major,minor,haplotypes)
	
	@classmethod
	def formatEntry(cls,entry):
		gt=entry.genotypes
		tof=[]
		tof.append(entry.chr)
		tof.append(str(entry.pos))
		tof.append(entry.refc)
		tof.append(entry.major+"/"+entry.minor)
		tof.append(" ".join(gt))
		form="\t".join(tof)
		return form
		

	