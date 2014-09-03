#!/usr/bin/env python


class PopulationHaplotype:
	"""
	Represents a line in a haplotype file, which corresponds to the genotypes of a whole population
	X       9121094 T       T/C     TC CT CC CT CT CC TT CT TT TC CC CC CC CT TT
	2L	300239	A	A/C	AA AA AA AA AA AA AA AA AA AA AA CC AA AA CC 
	"""
	def __init__(self,chrom,position,refchar,major,minor,haplotypes):	
		self.chrom=chrom
		self.position=position
		self.refchar=refchar
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
	def isPolymorphic(self):
		alleleset=set(self.haplotypes)
		assert(len(alleleset)>0)
		if(len(alleleset)>1):
			return True
		else:
			return False

class HaplotypeIO:
	
	@classmethod
	def genotypeCount(cls,inputFile):
		"""
		Return the number of genotypes present in a input file
		"""
		fh=open(inputFile)
		line=fh.next()
		fh.close()
		pa=HaplotypeIO.parseLine(line)
		return pa.genotypeCount
	
	@classmethod
	def haplotypeCount(cls,inputFile):
		fh=open(inputFile)
		line=fh.next()
		fh.close()
		pa=HaplotypeIO.parseLine(line)
		return pa.haplotypeCount
	

	@classmethod
	def parseLine(cls,line):
		a=line.split("\t");
		if(len(a)<5):
			raise StandardError("could not split line properly"+line)
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
		return PopulationHaplotype(a[0],a[1],a[2],major,minor,haplotypes)
	
	@classmethod
	def formatEntry(cls,entry):
		gt=entry.genotypes
		tof=[]
		tof.append(entry.chrom)
		tof.append(str(entry.position))
		tof.append(entry.refchar)
		tof.append(entry.major+"/"+entry.minor)
		tof.append(" ".join(gt))
		form="\t".join(tof)
		return form
		

	