import sys
import collections
from optparse import OptionParser, OptionGroup
import copy
import math
import random
from scipy.stats import poisson

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
	
	def subsample(self,targetcoverage):
		if targetcoverage ==0:
			return SinglePop(0,0,0,0,0,0)
		Af,Tf,Cf,Gf= (self.Af, self.Tf, self.Cf, self.Gf)
		assert(math.fabs(Af+Tf+Cf+Gf-1.0) <0.0000001)
		tuple=((Af,'A'),(Af+Tf,'T'),(Af+Tf+Cf,'C'),(1.0,'G'))
		
		Ac,Tc,Cc,Gc=(0,0,0,0)
		for c in range(0,targetcoverage):
			r=random.random()
			base=None
			for t in tuple:
				if r < t[0]:
					base=t[1]
					break
			if base=='A':
				Ac+=1
			elif base=='T':
				Tc+=1
			elif base=='C':
				Cc+=1
			elif base=='G':
				Gc+=1
			else:
				raise ValueError("Invalid subsampling")
		return SinglePop(Ac,Tc,Cc,Gc,0,0)

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
	def __init__(self,chr,pos,refc,populations):
		self.__chr=chr
		self.__pos=int(pos)
		self.__refc=refc
		self.__populations=populations

		
		

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
	
	@property
	def popcount(self):
		return len(self.__populations)
	
	def __str__(self):
		popstr="\t".join(map(str,self.populations))
		tojoin=[self.chr,self.pos,self.refc,popstr]
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
		
def get_mean_coverages(popcount,targetcoverage):
	toret=[]
	for i in range(0,popcount):
		pr=poisson.rvs(targetcoverage)
		toret.append(pr)
	return toret

def add_gc_bias(meancoverages,targetcoverage):
	rand=poisson.rvs(targetcoverage)
	cumprob=poisson.cdf(rand,targetcoverage) # cdf(x, mu, loc=0)	Cumulative density function.
	
	toret=[]
	for cov in meancoverages:
		if cov==0:
			toret.append(0)
		else:
			t=int(poisson.ppf(cumprob,cov)) # ppf(q, mu, loc=0)	Percent point function (inverse of cdf percentiles).
			toret.append(t)
	return toret
	
		
	


#Author: Dr. Robert Kofler
parser = OptionParser()
parser.add_option("--input", dest="input", help="the input file")
parser.add_option("--coverage",dest="coverage",help="the average coverage")
(options, args) = parser.parse_args()
overallcov=int(options.coverage)

meancoverages=None

for s in SyncReader(options.input):
	pops=s.populations
	if meancoverages is None:
		meancoverages=get_mean_coverages(len(pops),overallcov)
		
	poscoverage=add_gc_bias(meancoverages,overallcov)	
	
	nuevopops=[]
	for i,p in enumerate(pops):
		activecov=poscoverage[i]
		pr=0
		if(activecov>0):
			pr=poisson.rvs(activecov)
		np=p.subsample(pr)
		nuevopops.append(np)
	npl=PopLine(s.chr,s.pos,s.refc,nuevopops)
	print str(npl)

