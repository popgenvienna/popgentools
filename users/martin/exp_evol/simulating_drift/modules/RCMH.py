#!/usr/bin/env python
import sys,select
import collections
import CMH
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr


class Utility:
	
	@classmethod
	def getSinglePops(cls,pops,majora,minora):
		toret=[]
		for p in pops:
			(ac,tc,gc,cc)=(0,0,0,0)
			if(majora=="A"):
				ac=p[0]
			elif(majora=="T"):
				tc=p[0]
			elif(majora=="C"):
				cc=p[0]
			elif(majora=="G"):
				gc=p[0]
			else:
				raise StandardException("your fucked")
				
			if(minora=="A"):
				ac=p[1]
			elif(minora=="T"):
				tc=p[1]
			elif(minora=="C"):
				cc=p[1]
			elif(minora=="G"):
				gc=p[1]
			else:
				raise StandardException("your fucked")
			#(self,a,t,c,g,n,deletion)
			sp=CMH.SinglePop(ac,tc,cc,gc,0,0)
			toret.append(sp)
		return toret
			
	
	@classmethod
	def coverageValid(cls,populations,mincoverage,maxcoverage):
		if type(maxcoverage).__name__=='int':	
			for p in populations:
				if(p.cov < mincoverage or p.cov > maxcoverage):
					return False
			return True
		if type(maxcoverage).__name__=='list':
			for p in range(len(populations)):
				if(populations[p].cov < mincoverage or populations[p].cov > maxcoverage[p]):
					return False
			return True

	@classmethod
	def getDeletionCount(cls, populations):
		delcount=0
		for p in populations:
			delcount+=p.deletion
		return delcount
	
	@classmethod
	def getMajorAlleleCount(cls,populations):
		tal = {"A":0,"T":0,"C":0,"G":0}
		for p in populations:
			tal["A"]+=p.A
			tal["C"]+=p.C
			tal["T"]+=p.T
			tal["G"]+=p.G
			
		# transform into hash
		tar= [(key,value) for key,value in tal.items()]
		tar.sort(key=lambda x:-x[1])
		major=tar[0][0]
		minor=tar[1][0]
		
		toret=[]
		for p in populations:
			toret.append((p.countForAllele(major),p.countForAllele(minor)))
		return (toret,major,minor)
	
	@classmethod
	def numberAlleles(cls,populations):
		tal,sc = {"A":0,"T":0,"C":0,"G":0},0
		for p in populations:
			tal["A"]+=p.A
			tal["C"]+=p.C
			tal["T"]+=p.T
			tal["G"]+=p.G
			
		# transform into hash
		for k,v in tal.items():
			if v>0:
				sc+=1
		return sc		
		
		
	@classmethod
	def issnp(cls,populations,mincount):
		micount=0
		macount=0
		for p in populations:
			macount+=p[0]
			micount+=p[1]
		if micount>=mincount and macount>=mincount:
			return True
		else:
			return False
		
	@classmethod
	def getCMHPvalue(cls,populations):
		response=[]
		replicatecount=0
		for i in range(1,len(populations),2):
			base=populations[i-1]
			derived=populations[i]
			response.append(base[0])
			response.append(derived[0])
			response.append(base[1])
			response.append(derived[1])
			replicatecount+=1
		
		r=robjects.r
		response_robj = robjects.IntVector(response)
		dim_robj=robjects.IntVector([2,2,replicatecount])
		response_rar=robjects.r['array'](response_robj,dim=dim_robj)
		
		testres=r['mantelhaen.test'](response_rar)
		pvalue=testres[2][0]
		assert(pvalue<=1)
		assert(pvalue>=0)
		

		return pvalue

class SimulateDrift:
	
	@classmethod
	def multiSampleDrift(cls,popsize,alcount,derivedtimestamps,derivedcoverages):
		basecount=len(alcount)
		dtmatrix=SimulateDrift._createMatrix(derivedtimestamps,basecount)
		dcovmatrix=SimulateDrift._createMatrix(derivedcoverages,basecount)
		
		templist=[]
		for i in range(len(dtmatrix)):
			tbase=alcount[i]
			ttime=dtmatrix[i]
			tcov=dcovmatrix[i]
			tres=SimulateDrift.singleSampleDrift(popsize,tbase,ttime,tcov)
			templist.append(tres)
			
		toret=SimulateDrift._matrixToList(templist)
		return toret

	
	@classmethod
	def singleSampleDrift(cls,popsize,basepop,derivedtimestamps,derivedcoverages):
		#(1000,(50,50),(10,20,30,40),(30,40,60,70))
		if len(basepop)!=2:
			raise ValueError("Basepop has to have the size 2")
		if len(derivedtimestamps) != len(derivedcoverages):
			raise ValueError("length of the derivedtimes has to be equal to the derivedcoverages")
		if(sorted(derivedtimestamps) != sorted(derivedtimestamps)):
			raise ValueError("Derived time stamps have to be increasing")
		
		timedict=dict(zip(derivedtimestamps,derivedcoverages))
		
		majfreq=(float(basepop[0])/float(basepop[0]+basepop[1]) )
		ultimo_time=derivedtimestamps[-1]
		
		#create a set of the timestamps at which the tool has to report the allele frequencies
		stopset=set(derivedtimestamps)
		toret=[]
		for i in range(1, ultimo_time+1):
			majfreq=SimulateDrift._nextGeneration(popsize,majfreq)
			if i in stopset:
				cov=timedict[i]
				derfreq=SimulateDrift._sampleReads(cov,majfreq)
				toret.append(derfreq)
		return toret
	

	
	@classmethod
	def _sampleReads(cls,coverage,majorfreq):
		import random
		majcount=0
		mincount=0
		for i in range(coverage):
			r=random.random()
			if r < majorfreq:
				majcount+=1
			else:
				mincount+=1
		return (majcount,mincount)
	
	@classmethod
	def _nextGeneration(cls,popsize,majorfreq):
		import random
		majcount=0
		for i in range(popsize):
			r=random.random()
			# in case of 50%-50%  i want a balanced chance for every element
			# random is producing [0.0-1.0)
			# therefore [0.0-0.5) and [0.5-1.0) thus  r needs to be smaller than the major allele frequnecy
			if r < majorfreq:
				majcount+=1
		return float(majcount)/float(popsize)
	
	@classmethod
	def _createMatrix(cls,list,count):
		#(10,20,30,40,50,60),3 ->
		#((10,40),(20,50),(30,60))
		toret=[]
		for i in range(count):
			tlist=[]
			for k in range(i,len(list),3):
				tlist.append(list[k])
			toret.append(tlist)
		return toret
	
	@classmethod
	def _matrixToList(cls,matrix):
		#(10,40),(20,50),(30,60) ->
		#(10,20,30,40,50,60)
		
		toret=[]
		depth=len(matrix[0])
		for k in range(depth):
			for i in range(len(matrix)):
				toret.append(matrix[i][k])
		return toret
		