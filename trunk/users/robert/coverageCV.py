import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math


class SyncCovReader:
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
		covs=SyncCovReader.parse_sync_cov(a)
		return (chr,pos,covs)
	@classmethod
	def parse_sync_cov(cls,entries):
		eucov=[]
		for e in entries:
			a=map(float,e.split(":"))
			teuco=float(a[0])+float(a[1])+float(a[2])+float(a[3])
			eucov.append(teuco)
		return eucov

parser = OptionParser()

parser.add_option("--input",dest="input",help="the input file as sync")
(options, args) = parser.parse_args()

sH=collections.defaultdict(lambda:0)
ssH=collections.defaultdict(lambda:0)
nH=collections.defaultdict(lambda:0)
nunion=0;

for chr,pos,cov in SyncCovReader(options.input):
	coveredinany=False
	for i,c in enumerate(cov):
		if c>0:
			coverdinany=True
			sH[i]+=c
			ssH[i]+=c*c
			nH[i]+=1
	if coveredinany:
		nunion+=1


print "i\tmu\tstddev\tcv"
for k in sorted(sH.keys()):
	n=nH[k]
	mu=sH[k]/n
	muunion=sH[k]/nunion
	si= (ssH[k] - (sH[k]**2/n))/n
	s=math.sqrt(si)
	cv=s/mu
	print k,mu,muunion,si,cv,n




