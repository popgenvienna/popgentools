import sys
import collections



class SinglePop:
	def __init__(self,a,t,c,g,n,deletion):
		"""
		Represents the allele counts of a single population
		"""
		self.__a=int(a)
		self.__t=int(t)
		self.__c=int(c)
		self.__g=int(g)
		self.__n=int(n)
		self.__del=int(deletion)
		

	def issnp(self,mincount):
		"""
		>>> SinglePop(2,0,0,2,0,0).issnp(1)
		1
		>>> SinglePop(2,0,0,2,0,0).issnp(3)
		0
		>>> SinglePop(2,0,0,1,0,0).issnp(2)
		0
		>>> SinglePop(2,0,0,0,0,0).issnp(1)
		0
		>>> SinglePop(2,2,0,0,0,0).issnp(1)
		1
		"""
		alcount=self.count_alleles(mincount)
		if(alcount>1):
			return 1
		else:
			return 0
	

	def count_alleles(self,mincount):
		"""
		>>> SinglePop(2,0,0,2,0,0).count_alleles(1)
		2
		>>> SinglePop(2,0,0,2,0,0).count_alleles(2)
		2
		>>> SinglePop(2,0,0,2,0,0).count_alleles(3)
		0
		"""
		alcount=0
		if self.A>=mincount:
			alcount+=1
		if self.T>=mincount:
			alcount+=1
		if self.C>=mincount:
			alcount+=1
		if self.G>=mincount:
			alcount+=1
		return alcount
	
	def countForAllele(self,allele):
		"""
		>>> SinglePop(2,0,0,0,0,0).countForAllele("A")
		2
		>>> SinglePop(2,4,0,0,0,0).countForAllele("T")
		4
		"""
		return eval("self."+allele)		

	@property
	def A(self):
		"""
		>>> SinglePop(2,0,0,0,0,0).A
		2
		"""
		return self.__a
	
	@property
	def Af(self):
		"""
		>>> SinglePop(4,0,0,0,0,0).Af
		1.0
		>>> SinglePop(4,4,0,0,0,0).Af
		0.5
		"""
		return float(self.A)/float(self.cov)
	
	
	@property
	def T(self):
		"""
		>>> SinglePop(2,0,0,0,0,0).T
		0
		>>> SinglePop(0,2,0,0,0,0).T
		2
		"""
		return self.__t
	
	@property
	def Tf(self):
		"""
		>>> SinglePop(0,4,0,0,0,0).Tf
		1.0
		>>> SinglePop(4,4,0,0,0,0).Tf
		0.5
		"""
		return float(self.T)/float(self.cov)
		

	@property
	def C(self):
		"""
		>>> SinglePop(0,0,2,0,0,0).C
		2
		"""
		return self.__c
	
	@property
	def Cf(self):
		"""
		>>> SinglePop(0,0,4,0,0,0).Cf
		1.0
		>>> SinglePop(0,4,4,0,0,0).Cf
		0.5
		"""
		return float(self.C)/float(self.cov)

	@property
	def G(self):
		"""
		>>> SinglePop(0,0,0,3,0,0).G
		3
		"""
		return self.__g
	
	@property
	def Gf(self):
		"""
		>>> SinglePop(0,0,0,4,0,0).Gf
		1.0
		>>> SinglePop(0,4,0,4,0,0).Gf
		0.5
		"""
		return float(self.G)/float(self.cov)

	@property
	def N(self):
		"""
		>>> SinglePop(0,0,0,0,5,0).N
		5
		"""
		return self.__n
		
	
	@property
	def deletion(self):
		"""
		>>> SinglePop(0,0,0,0,0,6).deletion
		6
		"""
		return self.__del

	@property
	def cov(self):
		"""
		>>> SinglePop(2,3,1,1,0,0).cov
		7
		"""
		return self.A+self.T+self.C+self.G

	@property
	def totcov(self):
		"""
		>>> SinglePop(2,3,1,1,2,1).totcov
		10
		"""
		return self.cov+self.N+self.deletion
	
	def __str__(self):
		"""
		>>> str(SinglePop(6,5,4,3,2,1))
		'6:5:4:3:2:1'
		"""
		return ":".join(map(str,[self.A,self.T,self.C,self.G,self.N,self.deletion]))



class PopLine:
	def __init__(self,chr,pos,refc,populations,pvalue=None):
		self.__chr=chr
		self.__pos=int(pos)
		self.__refc=refc
		self.__populations=populations
		self.__pvalue=pvalue
		
		

	@property
	def chr(self):
		"""
		>>> PopLine("2L",1,"N",[],0.2).chr
		'2L'
		"""
		return self.__chr
	
	@property
	def pos(self):
		return self.__pos
	
	@property
	def refc(self):
		return self.__refc
	
	@property
	def populations(self):
		return self.__populations
	
	def subpopulations(self,populations):
		tpops=self.populations
		toret=[]
		for i in populations:
			toret.append(tpops[i-1])
		return toret
	
	@property
	def pvalue(self):
		return self.__pvalue
	
	@pvalue.setter
	def pvalue(self,value):
		self.__pvalue=value
	
	@property
	def popcount(self):
		return len(self.__populations)
	
	def __str__(self):
		"""
		"""
		popstr="\t".join(map(str,self.populations))
		pvalue=self.pvalue
		tojoin=[self.chr,self.pos,self.refc,popstr]
		if pvalue is not None:
			tojoin.append(pvalue)
		return "\t".join([str(x) for x in tojoin])
		

class SyncReader:
	def __init__(self,file):
		self.__filename=file
		if(isinstance(file,str)):
			self.__filehandle=open(file,"r")
		else:
			self.__filehandle=file

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
		population=[]
		for p in a:
			po=None
			if p=="-":
				po=SinglePop(0,0,0,0,0,0)
			else:
				s=p.split(":")
				po=SinglePop(*s)
			population.append(po)
		
		return PopLine(chr,pos,refc,population)

class CMHReader:
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
		pvalue=float(a.pop())
		
		population=[]
		for p in a:
			if p=="-":
				po=SinglePop(0,0,0,0,0,0)
			else:
				s=p.split(":")
				po=SinglePop(*s)
			population.append(po)
		
		return PopLine(chr,pos,refc,population,pvalue)
	




if __name__ == "__main__":
	import doctest
	doctest.testmod(verbose=1)


    

