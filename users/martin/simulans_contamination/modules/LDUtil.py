import sys
import collections
import numpy
import math
from CMH import CMHReader


	

class Util:
	
	#@classmethod
	#def meanstdv(x): ### calculate mean, median, stdev. standard error : x=datalist
	#	from math import sqrt
	#	n,mean,std,se,median = len(x),0,0,0,0
	#	for a in x:
	#	    mean = mean + a
	#	mean = mean / float(n)
	#	for a in x:
	#	    std = std + (a - mean)**2
	#	std = sqrt(std / float(n-1))
	#	se= std/sqrt(n)
	#	median=sorted(x)[len(x)/2]
	#	return mean,std,se,median
	
	@classmethod
	def geometricmean(cls,x):
		"""
		>>> Util.geometricmean([float("1e-12"),float("1e-14"),float("1e-10")])
		1.000000000000001e-12
		>>> Util.geometricmean([float("1e-2"),float("1e-2"),float("1e-14"),float("1e-10")])
		9.9999999999999943e-08
		"""
		if len(x)==0:
			return None
		logsum=sum([math.log(v) for v in x])
		logsum/=len(x)
		pval=math.exp(logsum)
		return pval
		
	@classmethod
	def mean(cls,x):
		"""
		>>> Util.geometricmean([float("1e-12"),float("1e-14"),float("1e-10")])
		1.000000000000001e-12
		>>> Util.geometricmean([float("1e-2"),float("1e-2"),float("1e-14"),float("1e-10")])
		9.9999999999999943e-08
		"""
		if len(x)==0:
			return None
		sum1=sum([v for v in x])
		pval=sum1/len(x)
		return pval	
		
	
	@classmethod
	def median(cls,x,log):
		"""
		>>> Util.median([float("1e-12"),float("1e-14"),float("1e-10")])
		9.9999999999999998e-13
		>>> Util.median([float("1e-11"),float("1e-20"),float("1e-12")])
		9.9999999999999998e-13
		>>> Util.median([float("1e-8"),float("1e-11"),float("1e-12"),float("1e-20")])
		3.1622776601683786e-12
		>>> Util.median([float("1e-8"),float("1e-10")])
		1.0000000000000007e-09
		"""
		mid =int(len(x)/2)
		sort=sorted(x)
		if len(x)==0:
			return None
		if len(x)%2==0:
			lower=sort[mid-1]
			upper=sort[mid]
			if log=="True":
				return Util.geometricmean([lower,upper])
			else:
				return Util.mean([lower,upper])
		else:
			return sort[mid]
	
	@classmethod
	def posh(x):
		h=collections.defaultdict(lambda:collections.defaultdict(lambda:0.0))
		for l in open(x,"r"):
			if l.rstrip()=="":
				continue
			p=PopIO.parse_cmhline(l)
			h[p.chr][p.pos]=p.pvalue
		return h

class CandidateSNP:
	def __init__(self,chr,pos,rangestart,rangeend,pvalue,ignorecandidate=1):
		self.__chr		= chr
		self.__pos		= pos
		self.__rangestart	= rangestart
		self.__rangeend		= rangeend
		self.__pvalue		= pvalue
		self.__neighborsnps	= []
		self.__ignorecandidate 	= ignorecandidate
		
	
	def appendSNP(self,snp):
		if(self.__ignorecandidate and self.chr == snp.chr and self.pos==snp.pos):
			return False
		self.__neighborsnps.append(snp)
		return True
		
	@property
	def chr(self):
		return self.__chr
	
	@property
	def pos(self):
		return self.__pos
		
	@property
	def pvalue(self):
		return self.__pvalue
	
	@property
	def snps(self):
		return self.__neighborsnps
	
	@property
	def rangeend(self):
		return self.__rangeend
	
	@property
	def rangestart(self):
		return self.__rangestart
	

	def distributeToBins(self,bincount):
		binar=[ [] for i in range(0,bincount)]
		dist=self.rangeend-self.rangestart+1
		binsize=int(dist/bincount)
		for snp in self.snps:
			arindex = snp.pos - self.rangestart
			arindex = int(arindex/binsize)
			binar[arindex].append(snp)
		return binar
			

class LDIO:
	
	@classmethod
	def read_candidatehash(cls,filename,maxrange):
		
		chrh=collections.defaultdict(lambda:collections.defaultdict(lambda:[]))
		candl=[]
		for p in CMHReader(filename):
			chr=p.chr
			pos=p.pos
			rangestart=pos-maxrange
			rangend=pos+maxrange-1
			cand=CandidateSNP(p.chr,p.pos,rangestart,rangend,p.pvalue)
			candl.append(cand)
			for i in range(rangestart,rangend+1):
				chrh[chr][i].append(cand)
		return chrh,candl


if __name__ == "__main__":
	import doctest
	doctest.testmod(verbose=1)


    

