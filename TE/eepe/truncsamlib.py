import re


class PTruncCollection:
	def __init__(self,sams):
		self.sams=sams
	
	def sams(self):
		return self.sams
	
	def get_truncations(self):
		truncs=[]
		for s in self.sams:
			truncs.extend(s.get_truncations())
		return truncs
	
	def get_genomicfragments(self):
		frags=[]
		for s in self.sams:
			frags.extend(s.get_genomefragments())
		return frags
	
	def get_coverage(self):
		tmp=self.get_genomicfragments()
		tmp.extend(self.get_truncations())
		cov=[0]*2908
		for i in tmp:
			start=i[0]
			end=i[1]
			for k in range(start,end+1):
				cov[k]+=1
		return cov
	
	def get_averagecoverage(self,bound):
		cov=self.get_coverage()
		start=bound
		end=len(cov)-bound
		cs=0
		c=0
		for i in range(start,end):
			c+=1
			cs+=cov[i]
		return float(cs)/float(c)




class PTruncSamEntry:
	@classmethod
	def load(cls,file):
		toret=[]
		for l in open(file):
			if l.startswith("@"):
				continue
			a=l.rstrip("\n").split("\t")
			# FCC4M7EACXX:2:2101:4216:63007#CGATGTAT	0	PPI251	1	40	100M	*	0	0	CATGATGAAATAACATAAGGTGGTCCCGTCGAAAGCCGAAGCTTACCGAAGTATACACTTAAATTCAGTGCACGTTTGCTTGTTGAGAGGAAAGGTTGTG	@@@DDDDDDBHFGIDGIIIG<CD@EGHIFGHGEHIGGIIEGHHCDGHG8B;CAHHAHEHFFFFFFEDECCCCDDDBDDDDDDDDD@A<C?8<9AB(:?B<	MD:Z:100	NH:i:1	HI:i:1	NM:i:0	SM:i:40	XQ:i:40	X2:i:0	XO:Z:UU	XG:Z:A
			
			chr=a[2]
			start=int(a[3])
			cig=a[5]
			if cig=="*" or chr=="*":
				continue
			s=PTruncSamEntry(chr,start,cig)
			toret.append(s)
			
		return PTruncCollection(toret)	
	
	
	def __init__(self,chr,start,cigar):
		if chr !="PPI251":
			raise Exception("read not mapping to the P-element {0}".format(chr))
		self.__start=start
		self.__cigar=cigar
		pe=[]
		for fi in re.finditer(r"(\d+)([HSIDMN])", cigar):
			num=int(fi.group(1))
			id=fi.group(2)
			pe.append((num,id))
		self.__splitcigar=pe	
		
	def cigar(self):
		return self.__cigar
		
	def start(self):
		return self.__start


	def get_truncations(self):
		end=self.__start
		cigsplit=self.__splitcigar
		truncations=[]
		for num,id in cigsplit:
			if id=="M":
				end+=num
			elif id=="D":
				end+=num
			elif id=="N":
				s=end
				end+=num
				e=end-1
				truncations.append((s,e))
			elif id=="I" or id=="S" or id=="H":
				pass
			else:
				raise Exception("unknown cigar"+id)
		return truncations

	def get_genomefragments(self):
		end=self.__start
		cigsplit=self.__splitcigar
		fragments=[]
		for num,id in cigsplit:
			if id=="M":
				s=end
				end+=num
				e=end-1
				fragments.append((s,e))
			elif id=="D":
				s=end
				end+=num
				e=end-1
				fragments.append((s,e))
			elif id=="N":
				end+=num
			elif id=="I" or id=="S" or id=="H":
				pass
			else:
				raise Exception("unknown cigar"+id)
		return fragments

	def end(self):
		end=self.__start
		cigsplit=self.__splitcigar
		for num,id in cigsplit:
			if id=="M":
				end+=num
			elif id=="D":
				end+=num
			elif id=="N":
				end+=num
			elif id=="I" or id=="S" or id=="H":
				pass
			else:
				raise Exception("unknown cigar"+id)
		return end-1

class PHandling:
	

	
	@classmethod
	def getavcov_bounded(cls,cov,start,end):
		cs=0
		c=0
		for i in range(start,end+1):
			c+=1
			cs+=cov[i]
		return float(cs)/float(c)


