import sys
import collections
from optparse import OptionParser, OptionGroup
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python allele_freq_changes.py -i candidates.cmh -b 100 -p 1,3,5 -o out
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	The input file (-c) can be either a *.sync or a cmh test output. This script calculates the allele frequency changes of the minor frequent allele based on the first population set in the list of populations to compare (-p), e.g. population 1 if parameter is set as -p 1,3,5. If this population is fixed for one allele (i.e. min. allele freq: 1.0) the other populations will be searched for the nucleotide of the minor allele frequency. If all three populations are fixed for one allele, allele frequency changes will be 0 for all combinations. Then, the changes in frequencies of this nucleotide will be estimated in a pairwise manner for all population combinations by substracting the allele-frequencies of the second (in our case: the later population) in the pairwise comparison from the value of the first. population. The script produces four different outputs: 1) out.freqchanges: a table with all frequency changes for every position separately. The values can be either positive=increase in frequency, or negative=decrease in frequency. 2) out.summary: a table for all populations summarizing the absolute changes for all positions of the input file (minimum value, mean, median, 1st and 3rd quantile, maximum) 3) out.boxplot: a Boxplot showing the absolute frequency changes of the individual comparisons. The whiskers reach until the 1st or 3rd quantile +/- 1.5 x interquantile range (see R documentation). 4) out.histXXX.pdf: Separte histograms for all pairwise comparison showing the distribution of allele frequency changes.  
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-b", "--bins", dest="bins", help="No. of bins for the histograms")
parser.add_option("-o", "--out", dest="out", help="outputfile containing SNPs, which show high freq changes starting from very low (0) frequencies")
parser.add_option("-p", "--pops", dest="pops", help="define 2 or 3 populations for which you want to calculate the allele frequency change. populations are separated by a \",\"(e.g. 1,2)")
parser.add_option("-n", "--pnames", dest="pnames", help="define names of 2 or 3 populations for which you want to calculate the allele frequency change. populations are separated by a \",\"(e.g. 1,2)")
parser.add_option_group(group)
(options, args) = parser.parse_args()

#2L	4910	A	6:0:0:0:0:0	-	5:0:0:0:0:0	25:0:0:0:0:0	33:0:0:0:0:0
#2L	4911	G	0:0:0:7:0:0	-	0:0:0:5:0:0	0:0:0:27:0:0	0:0:0:33:0:0
#2L	4912	A	7:0:0:0:0:0	-	5:0:0:0:0:0	27:0:0:0:0:0	32:0:0:0:0:0
#2L	4913	G	0:0:0:7:0:0	-	0:0:0:5:0:0	0:0:0:27:0:0	0:0:0:34:0:0
#2L	4914	A	7:0:0:0:0:0	-	5:0:0:0:0:0	26:0:0:0:0:0	33:0:0:0:0:0
#2L	4915	G	0:0:0:8:0:0	-	0:0:0:5:0:0	0:0:0:22:0:0	0:0:0:33:0:0
#2L	4916	A	8:0:0:0:0:0	-	6:0:0:0:0:0	25:0:0:0:0:0	33:0:0:0:0:0
#2L	4917	G	0:0:0:8:0:0	-	0:0:0:6:0:0	0:0:0:23:0:0	0:0:0:32:0:0
#2L	4918	C	0:0:8:0:0:0	-	0:0:6:0:0:0	0:0:22:0:0:0	0:0:30:0:0:0

def allelefreqhash(inp,pops):
	allelefreqchange={}
	sc=0
	ha,empty=[],[]
	a1,b1=str(pops).split(","),str(pops).split(",")
	for i in a1:
		for j in b1[1:]:
			if j!=i:
				ha.append(str(i)+","+str(j))
				empty.append(0)	
		a1=a1[1:]	
		b1=b1[1:]
	sync=open(str(inp),"r")
	for l in sync:
		a=l.split()
		if "-" not in a[:-1]:
			b=a[0]+"_"+a[1]
			allelefreqchange[b]= {}.fromkeys(ha,0) #dict(zip(ha,empty))
	sync=open(str(inp),"r")
	for l in sync:
		fr2,fr2a={},{}
		sc+=1
		a=l.split()
		b=a[0]+"_"+a[1]
		if ":" not in l.split()[-1:]:
			if "-" not in a[:-1] and "\t0:0:0:0:0:0" not in a[:-1]:
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
		if ":" in l.split()[-1:]:
			if "-" not in l and "\t0:0:0:0:0:0" not in l:
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
					
	return allelefreqchange
	
	
allelefreqchange=allelefreqhash(str(options.inp),str(options.pops))
								
li,key=[],[]
allelefreqlistabs,allelefreqlist=collections.defaultdict(lambda:[]),collections.defaultdict(lambda:[])
c=0
out1=open(str(options.out)+".freqchanges","w")
for k,v in allelefreqchange.items():
	for k1,v1 in v.items():
		key.append(k1)
	if c==0:
		out1.write("chrom\tstart\t"+"\t".join(key)+"\n")
		c=1
	if c==1:
		for k1,v1 in v.items():
			li.append(v1)
			allelefreqlistabs[k1].append(abs(float(v1)))
			allelefreqlist[k1].append(float(v1))
		out1.write(k.split("_")[0]+"\t"+k.split("_")[1]+"\t"+"\t".join(map(str,li))+"\n")
		li=[]
popnames=dict(zip(str(options.pops).split(","),str(options.pnames).split(",")))
if len(str(options.pops).split(","))==4:
	li,key,li1,key1=[],[],[],[]
	r=robjects.r
	graphics = importr('graphics')
	grdevices = importr('grDevices')
	
	for k,v in sorted(allelefreqlistabs.items()):
		key.append(popnames[k.split(",")[0]]+","+popnames[k.split(",")[1]])
		li.append(v)
	for k2,v2 in sorted(allelefreqlist.items()):
		li1.append(v2)
	names = robjects.StrVector(key)
	title=robjects.StrVector(str(options.inp).split("/")[-1:])
	values1= robjects.vectors.FloatVector(li[0])
	values2= robjects.vectors.FloatVector(li[1])
	values3= robjects.vectors.FloatVector(li[2])
	values5= robjects.vectors.FloatVector(li[3])
	values6= robjects.vectors.FloatVector(li[4])
	values4= robjects.vectors.FloatVector(li[5])
	values1a= robjects.vectors.FloatVector(li1[0])
	values2a= robjects.vectors.FloatVector(li1[1])
	values3a= robjects.vectors.FloatVector(li1[2])
	values4a= robjects.vectors.FloatVector(li1[3])
	values5a= robjects.vectors.FloatVector(li1[4])
	values6a= robjects.vectors.FloatVector(li1[5])
	
	grdevices.pdf(file=str(options.out)+'_boxplot.pdf',width=20,height=10)
	graphics.boxplot(values1,values2,values3,values4,values5,values6,names=names,main=title,ylim=r.range(0,1))
	grdevices.dev_off()
	out2=open(str(options.out)+".summary","w")
	out2.write( "Comparison\tMin.\t1st_Qu.\tMedian\tMean\t3rd_Qu.\tMax\n")
	out2.write(names[0]+"\t"+str(r.summary(values1)[0])+"\t"+str(r.summary(values1)[1])+"\t"+str(r.summary(values1)[2])+"\t"+str(r.summary(values1)[3])+"\t"+str(r.summary(values1)[4])+"\t"+str(r.summary(values1)[5])+"\n")
	out2.write(names[1]+"\t"+str(r.summary(values2)[0])+"\t"+str(r.summary(values2)[1])+"\t"+str(r.summary(values2)[2])+"\t"+str(r.summary(values2)[3])+"\t"+str(r.summary(values2)[4])+"\t"+str(r.summary(values2)[5])+"\n")
	out2.write(names[2]+"\t"+str(r.summary(values3)[0])+"\t"+str(r.summary(values3)[1])+"\t"+str(r.summary(values3)[2])+"\t"+str(r.summary(values3)[3])+"\t"+str(r.summary(values3)[4])+"\t"+str(r.summary(values3)[5])+"\n")
	out2.write(names[3]+"\t"+str(r.summary(values4)[0])+"\t"+str(r.summary(values4)[1])+"\t"+str(r.summary(values4)[2])+"\t"+str(r.summary(values4)[3])+"\t"+str(r.summary(values4)[4])+"\t"+str(r.summary(values4)[5])+"\n")
	out2.write(names[4]+"\t"+str(r.summary(values5)[0])+"\t"+str(r.summary(values5)[1])+"\t"+str(r.summary(values5)[2])+"\t"+str(r.summary(values5)[3])+"\t"+str(r.summary(values5)[4])+"\t"+str(r.summary(values5)[5])+"\n")
	out2.write(names[5]+"\t"+str(r.summary(values6)[0])+"\t"+str(r.summary(values6)[1])+"\t"+str(r.summary(values6)[2])+"\t"+str(r.summary(values6)[3])+"\t"+str(r.summary(values6)[4])+"\t"+str(r.summary(values6)[5])+"\n")
	
	grdevices.pdf(str(options.out)+'_hist'+names[0]+'.pdf')
	graphics.hist(values1a,breaks=int(options.bins),main=names[0],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
	grdevices.pdf(str(options.out)+'_hist'+names[1]+'.pdf')
	graphics.hist(values2a,breaks=int(options.bins),main=names[1],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
	grdevices.pdf(str(options.out)+'_hist'+names[2]+'.pdf')
	graphics.hist(values3a,breaks=int(options.bins),main=names[2],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
	grdevices.pdf(str(options.out)+'_hist'+names[3]+'.pdf')
	graphics.hist(values4a,breaks=int(options.bins),main=names[3],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
	grdevices.pdf(str(options.out)+'_hist'+names[4]+'.pdf')
	graphics.hist(values5a,breaks=int(options.bins),main=names[4],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
	grdevices.pdf(str(options.out)+'_hist'+names[5]+'.pdf')
	graphics.hist(values6a,breaks=int(options.bins),main=names[5],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
if len(str(options.pops).split(","))==3:
	li,key,li1,key1=[],[],[],[]
	r=robjects.r
	graphics = importr('graphics')
	grdevices = importr('grDevices')
	
	for k,v in sorted(allelefreqlistabs.items()):
		key.append(popnames[k.split(",")[0]]+","+popnames[k.split(",")[1]])
		li.append(v)
	for k2,v2 in sorted(allelefreqlist.items()):
		li1.append(v2)
	names = robjects.StrVector(key)
	title=robjects.StrVector(str(options.inp).split("/")[-1:])
	values1= robjects.vectors.FloatVector(li[0])
	values2= robjects.vectors.FloatVector(li[1])
	values3= robjects.vectors.FloatVector(li[2])
	values1a= robjects.vectors.FloatVector(li1[0])
	values2a= robjects.vectors.FloatVector(li1[1])
	values3a= robjects.vectors.FloatVector(li1[2])
	
	grdevices.pdf(str(options.out)+'_boxplot.pdf')
	graphics.boxplot(values1,values2,values3,names=names,main=title,ylim=r.range(0,1))
	grdevices.dev_off()
	out2=open(str(options.out)+".summary","w")
	out2.write( "Comparison\tMin.\t1st_Qu.\tMedian\tMean\t3rd_Qu.\tMax\n")
	out2.write(names[0]+"\t"+str(r.summary(values1)[0])+"\t"+str(r.summary(values1)[1])+"\t"+str(r.summary(values1)[2])+"\t"+str(r.summary(values1)[3])+"\t"+str(r.summary(values1)[4])+"\t"+str(r.summary(values1)[5])+"\n")
	out2.write(names[1]+"\t"+str(r.summary(values2)[0])+"\t"+str(r.summary(values2)[1])+"\t"+str(r.summary(values2)[2])+"\t"+str(r.summary(values2)[3])+"\t"+str(r.summary(values2)[4])+"\t"+str(r.summary(values2)[5])+"\n")
	out2.write(names[2]+"\t"+str(r.summary(values3)[0])+"\t"+str(r.summary(values3)[1])+"\t"+str(r.summary(values3)[2])+"\t"+str(r.summary(values3)[3])+"\t"+str(r.summary(values3)[4])+"\t"+str(r.summary(values3)[5])+"\n")
	
	grdevices.pdf(str(options.out)+'_hist'+names[0]+'.pdf')
	graphics.hist(values1a,breaks=int(options.bins),main=names[0],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
	grdevices.pdf(str(options.out)+'_hist'+names[1]+'.pdf')
	graphics.hist(values2a,breaks=int(options.bins),main=names[1],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
	grdevices.pdf(str(options.out)+'_hist'+names[2]+'.pdf')
	graphics.hist(values3a,breaks=int(options.bins),main=names[2],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
if len(str(options.pops).split(","))==2:
	li,key,li1,key1=[],[],[],[]
	r=robjects.r
	graphics = importr('graphics')
	grdevices = importr('grDevices')
	for k,v in sorted(allelefreqlistabs.items()):
		key.append(popnames[k.split(",")[0]]+","+popnames[k.split(",")[1]])
		li.append(v)
	for k2,v2 in sorted(allelefreqlist.items()):
		li1.append(v2)
	names = robjects.StrVector(key)
	values1= robjects.vectors.FloatVector(li[0])
	values2= robjects.vectors.FloatVector(li[1])
	values1a= robjects.vectors.FloatVector(li1[0])
	values2a= robjects.vectors.FloatVector(li1[1])
	
	grdevices.pdf(str(options.out)+'_boxplot.pdf')
	graphics.boxplot(values1,values2,names=names,ylim=r.range(0,1))
	grdevices.dev_off()
	out2=open(str(options.out)+".summary","w")
	out2.write( "Comparison\tMin.\t1st_Qu.\tMedian\tMean\t3rd_Qu.\tMax\n")
	out2.write(names[0]+"\t"+str(r.summary(values1)[0])+"\t"+str(r.summary(values1)[1])+"\t"+str(r.summary(values1)[2])+"\t"+str(r.summary(values1)[3])+"\t"+str(r.summary(values1)[4])+"\t"+str(r.summary(values1)[5])+"\n")
	out2.write(names[1]+"\t"+str(r.summary(values2)[0])+"\t"+str(r.summary(values2)[1])+"\t"+str(r.summary(values2)[2])+"\t"+str(r.summary(values2)[3])+"\t"+str(r.summary(values2)[4])+"\t"+str(r.summary(values2)[5])+"\n")
	grdevices.pdf(str(options.out)+'_hist'+names[0]+'.pdf')
	graphics.hist(values1a,breaks=int(options.bins),main=names[0],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
	grdevices.pdf(str(options.out)+'_hist'+names[1]+'.pdf')
	graphics.hist(values2a,breaks=int(options.bins),main=names[1],xlab="",xlim=r.range(-0.5,0.5))
	grdevices.dev_off()
