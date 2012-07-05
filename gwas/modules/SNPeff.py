import sys
import collections
import re



class SNPeff:
	def __init__(self,chr,pos,refc,change,geneid,genename,transcriptid,effect,line):
		self.__chr=chr
		self.__pos=pos
		self.__refc=refc
		self.__change=change
		self.__geneid=geneid
		self.__genename=genename
		self.__transcriptid=transcriptid
		self.__effect=effect
		self.__line=line
	
	def __str__(self):
		return self.__line
	
	@property
	def geneid(self):
		return self.__geneid
	
	@property
	def genename(self):
		return self.__genename

	@property
	def chr(self):
		return self.__chr
	
	@property
	def pos(self):
		return self.__pos
	
	@property
	def effect(self):
		return self.__effect

class SNPeffReader:
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")
	
	def __iter__(self):
		return self
	
	def next(self):
		line=self.__filehandle.readline()
		if re.match("# Chrom",line):
			line=self.__filehandle.readline()
		if line=="":
			raise StopIteration
		line=line.rstrip('\n')
		a=line.split("\t")
		chrom=a[0]
		pos=int(a[1])
		refc=a[2]
		change=a[3]
		geneid=a[9]
		genename=a[10]
		transcriptid=a[12]
		fulleffect=a[15]
		effect=re.search(r"^([^:]+)",fulleffect).group(1)
		return SNPeff(chrom,pos,refc,change,geneid,genename,transcriptid,effect,line)


if __name__ == "__main__":
	import doctest
	doctest.testmod(verbose=1)


    

