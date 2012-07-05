import collections
import sys 
from optparse import OptionParser, OptionGroup
import time
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import math

#Author: Martin Kapun
#version 1.0

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python stauffer.py --inp data1.fet+data2.fet --pops 1,2+1,4 --sync fulldata.sync --out data.stauffer  
2)	the rpy2 package needs to be installed: http://sourceforge.net/projects/rpy/files/rpy2/ and type \"sudo python setup.py install\" in the directory
	""")

parser.add_option("-o", "--out", dest="out",help="output file")
parser.add_option("-i", "--inp", dest="inp",help="output from Ram's FET script. The files must be separated by a \"+\"")
parser.add_option("-p", "--pops", dest="pops",help="ID of the population in the FET comparison. Must be linked to the position of the corresponding columns in the *.sync file. E.g. \"1,2\" stands for the FET pvalues of the first two populations in the sync file. The ID's must also be separated by a \"+\" (e.g. 1,2+1,4)")
parser.add_option("-s", "--sync", dest="sync", help="synchronized file")
parser.add_option_group(group)
(options, args) = parser.parse_args()	
	
#########################################################   CODE   #########################################################################

inp=(str(options.inp).split("+"))
pop=(str(options.pops).split("+"))
fulllist=[]
fetdict={}
c=0
for l in open(inp[0],"r"):
	a=l.split()
	b=a[0]+"_"+a[1]
	fetdict[b]= {}.fromkeys(pop,0.0)
	c+=1
	if c%1000000==0:
		print str(c)+" lines filled in hash: fetdict"

for i in range(0,len(inp)):
	list1=[]
	ids=[]
	c=0
	for l in open(inp[i],"r"):
		a=l.split()
		b=a[0]+"_"+a[1]
		if b in fetdict:
			fetdict[b][pop[i]]=float(a[2])
			c+=1
		if c%100000==0:
			print str(c)+" lines processed in: "+str(inp[i])
#print fetdict
fet2dict={}
c=0
count=0
for k,v in fetdict.items():
	#print v
	if sum(v.values())!=len(inp) and 0.0 not in v.values() :
		fet2dict[k]=v
		count+=1
		if count%100000==0:
			print str(count)+" lines filled in new hash: fet2dict"
fetdict={}
#print 
out=open(str(options.out),"w")	
out.write("""Chromosome	position	"""+"\t".join(pop)+"""	Stauffer-score	product_of_pvalues"""+"\n")
for l in open(str(options.sync),"r"):
		fr2,fr2a={},{}
		a=l.split()
		b=a[0]+"_"+a[1]
		if "-" not in a and "\t0:0:0:0:0:0" not in a[:-1]:
			pop=a[int(str(pops).split(",")[0])+2]
			bases=["A","T","C","G"]
			countdict=dict(zip(bases,map(float,pop.split(":")[:4])))
			coverage=sum(countdict.values())
			frequencies=(str(float(countdict["A"])/coverage)+":"+str(float(countdict["T"])/coverage)+":"+str(float(countdict["C"])/coverage)+":"+str(float(countdict["G"])/coverage))
			freqdict=dict(zip(bases,map(float,frequencies.split(":")[:4])))
			for k,v in freqdict.items():
				if v!=0:
					fr2[k]=v
			minallele=min(fr2,key=lambda a: fr2.get(a))
			freqch=0
			if freqdict[minallele]!=1.0:
				a1=str(pops).split(",")
				b1=str(pops).split(",")
				for i in a1:
					for j in b1[1:]:
						if j!=i:
							comp1=a[int(i)+2]
							comp2=a[int(j)+2]
							f1=dict(zip(bases,map(float,comp1.split(":")[:4])))
							f2=dict(zip(bases,map(float,comp2.split(":")[:4])))
							cov1=sum(f1.values())
							cov2=sum(f2.values())
							if cov1!=0 and cov2!=0:
								freq1=(str(float(f1["A"])/cov1)+":"+str(float(f1["T"])/cov1)+":"+str(float(f1["C"])/cov1)+":"+str(float(f1["G"])/cov1))
								fr1=dict(zip(bases,map(float,freq1.split(":")[:4])))
								freq2=(str(float(f2["A"])/cov2)+":"+str(float(f2["T"])/cov2)+":"+str(float(f2["C"])/cov2)+":"+str(float(f2["G"])/cov2))
								fr2=dict(zip(bases,map(float,freq2.split(":")[:4])))
								allelefreqchange[b][str(i)+","+str(j)]=(float(fr2[minallele])-float(fr1[minallele]))
							else:
								continue
							
					a1=a1[1:]	
					b1=b1[1:]
			if freqdict[minallele]==1.0:
				minkey2="A"			
				minval=1.0
				for i in str(pops).split(","):
					fr2a={}
					if minval==1.0:
						comp=a[int(i)+2]
						f=dict(zip(bases,map(float,comp.split(":")[:4])))
						cov=sum(f.values())
						if cov!=0:
							freq=(str(float(f["A"])/cov)+":"+str(float(f["T"])/cov)+":"+str(float(f["C"])/cov)+":"+str(float(f["G"])/cov))
							fr=dict(zip(bases,map(float,freq.split(":")[:4])))
							for k,v in fr.items():
								if v!=0:
									fr2a[k]=v
							minkey2=min(fr2a,key=lambda a: fr2a.get(a))
							minval=float(fr[minkey2])
						else:
							continue
					else:
						break
						minallele=minkey2
				minallele=minkey2
		if b in fet2dict:
			#print fet2dict[b].keys()
			c+=1
			if c%100000==0:
				print str(c)+" lines processed of list of FET p-values<1.0 (count:"+str(len(fet2dict))+")"
			r=robjects.r
			subfunction=0
			result=0
			pvalueprod=1
			pvlist=[]
			for i in range(0,len(fet2dict[b])):
				pvalues=fet2dict[b].values()
				popindices=fet2dict[b].keys()
				synccov=[]
				for line in l.split():
					if ":" in line:
						synccov.append(sum(map(int, line.split(":")[:4])))
				#print synccov,popindices
				ni=float(synccov[int(popindices[i].split(",")[0])-1])/(float(synccov[int(popindices[i].split(",")[0])-1])+float(synccov[int(popindices[i].split(",")[1])-1]))
				wi=math.sqrt(ni)
				if pvalues[i]=="1":
					pv=float(0.9999999999999999)
				else:
					pv=float(pvalues[i])
				subfunction+=float(wi)*float(r['qnorm'](float(1-float(pv)))[0])
				pvalueprod*=pv
				pvlist.append(str(pv))
			result=1-r['pnorm'](subfunction)[0]
			out.write("\t".join(b.split("_"))+"\t"+"\t".join(pvlist)+"\t"+str(result)+"\t"+str(pvalueprod)+"\n")
#	if "-" in l:
#		a=l.split()
#		print a[0]+"\t"+a[1]+"\tstauffer score not available, too less coverage"