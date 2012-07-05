import sys
import collections
from optparse import OptionParser, OptionGroup
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.functions import SignatureTranslatedFunction
from optparse import OptionParser,OptionGroup
import copy

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python allele_freq_boxplot.py - candidates_BF15.cmh,candidates_BF37.cmh,candidates_F15F37.cmh -p 1,3,5 -c BF15,BF37,F15F37 -o out
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script produces boxplots of three datasets based on our cmh comparison  
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-s", "--sync", dest="sync", help="*.sync or cmh output file")
parser.add_option("-p", "--pops", dest="pops", help="define 2 or 3 populations for which you want to calculate the allele frequency change. populations are separated by a \",\"(e.g. 1,2)")
parser.add_option("-n", "--pnames", dest="pnames", help="define names of 2 or 3 populations for which you want to calculate the allele frequency change. populations are separated by a \",\"(e.g. 1,2)")
parser.add_option("-c", "--cand", dest="cand", help="define names of 2 or 3 populations for which you want to calculate the allele frequency change. populations are separated by a \",\"(e.g. 1,2)")
parser.add_option("-o", "--out", dest="out", help="define names of 2 or 3 populations for which you want to calculate the allele frequency change. populations are separated by a \",\"(e.g. 1,2)")
parser.add_option_group(group)
(options, args) = parser.parse_args()

def combinations(pops): 
	ha={}
	a1,b1=pops,pops
	for i in a1:
		for j in b1[1:]:
			if j!=i:
				ha[(str(i)+","+str(j))]=[]
		a1=a1[1:]	
		b1=b1[1:]
	return ha

def minallele(inp,pops):
	h,h2=collections.defaultdict(lambda:0),{}
	for i in range(0,len(pops)):
		h["A"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[0])
		h["T"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[1])
		h["C"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[2])
		h["G"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[3])
	cov=sum(h.values())
	#print l, h,cov
	if cov!=0:
		for k,v in h.items():
			if float(v)/cov!=0.0:
				h2[k]=float(v)/cov
		minallele=min(h2,key=lambda a: h2.get(a))
		#print max(h2.values())
		if max(h2.values())!=1.0:
			return minallele
		else:
			return "na"
	else:
		return "na"
		
def meanstdv(x):
    from math import sqrt
    n, mean, std,se = len(x), 0, 0,0
    for a in x:
		mean = mean + a
    mean = mean / float(n)
    for a in x:
		std = std + (a - mean)**2
    std = sqrt(std / float(n-1))
    se=std/sqrt(n)
    return mean, std,se

	
def freqhash(inp,i):
	h,h2=collections.defaultdict(lambda:0),{}
	h["A"]+=int(inp.split()[i+2].split(":")[0])
	h["T"]+=int(inp.split()[i+2].split(":")[1])
	h["C"]+=int(inp.split()[i+2].split(":")[2])
	h["G"]+=int(inp.split()[i+2].split(":")[3])
	cov=sum(h.values())
	if cov!=0:
		for k,v in h.items():
			h2[k]=float(v)/cov
		#print h2
		return h2
	else:
		return "na"
		
def max_allele_change(inp,minpop,maxpop):
	h=collections.defaultdict(lambda:0)
	pop1=freqhash(inp,int(minpop))
	pop2=freqhash(inp,int(maxpop))
	for k,v in pop1.items():
		h[k]=pop2[k]-v
	maxallele=max(h,key=lambda a: h.get(a))
	#	if h[maxallele]<0.05:
	#		print inp.rstrip()+"\t"+str(maxallele)+"\t"+str(minpop)+"\t"+str(maxpop)+"\t"+str(h[maxallele])
	return maxallele
out=str(options.out)
pop=str(options.pnames).split(",")
if len(str(options.pops).split(","))==4:
	pops=str(options.pops).split(",")
	cand=str(options.cand).split(",")
	pop=str(options.pnames).split(",")
	fullhash,newhash={},{}
	feat=["mean","stdabw","se"]
	for l in cand:
		fullhash[l]=dict(zip(pops,[[],[],[],[]]))
		newhash[l]=dict(zip(feat,[[],[],[]]))
	for l in open(str(options.sync).split(",")[0],"r"):
		for i in range(0,len(pops)):
			#print l.rstrip()+"\t"+str(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))])+"\t"+str(max_allele_change(l,int(pops[0]),int(pops[1])))+"\t"+str(pops[1])
			fullhash["bf15"][pops[i]].append(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))])
			#print freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))]
	for l in open(str(options.sync).split(",")[1],"r"):
		for i in range(0,len(pops)):
			#print l.rstrip()+"\t"+str(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))])+"\t"+str(max_allele_change(l,int(pops[0]),int(pops[1])))+"\t"+str(pops[1])
			fullhash["bf37"][pops[i]].append(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[3]))])
			#print fullhash
	for l in open(str(options.sync).split(",")[2],"r"):
		for i in range(0,len(pops)):
			#print l.rstrip()+"\t"+str(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))])+"\t"+str(max_allele_change(l,int(pops[0]),int(pops[1])))+"\t"+str(pops[1])
			fullhash["f15f37"][pops[i]].append(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[1]),int(pops[3]))])
	#print fullhash
	for k,v in sorted(fullhash.items()):
		for k1,v1 in sorted(v.items()):
			newhash[k]["mean"].append(meanstdv(v1)[0])
	for k,v in sorted(fullhash.items()):
		for k1,v1 in sorted(v.items()):
			newhash[k]["stdabw"].append(meanstdv(v1)[1])
	for k,v in sorted(fullhash.items()):
		for k1,v1 in sorted(v.items()):
			newhash[k]["se"].append(meanstdv(v1)[2])
	#print newhash
	
	r=robjects.r
	graphics = importr('graphics')
	grdevices = importr('grDevices')
	gplot=importr('gplots')
	graphics.par = SignatureTranslatedFunction(graphics.par,init_prm_translate = {'cex_axis': 'cex.axis'})
	graphics.par1 = SignatureTranslatedFunction(graphics.par,init_prm_translate = {'cex_lab': 'cex.lab'})
	popnames = robjects.StrVector(pop)
	print popnames
	#sys.exit()
	grdevices.pdf(file=out,width=30,height=10)
	a_bf15=[0.5,1.67,2.84,4.01]
	candr=robjects.vectors.StrVector(cand)
	colr=robjects.vectors.StrVector(["green","red","blue"])
	
	base_bf15=robjects.vectors.FloatVector(fullhash["bf15"][pops[0]])
	f15_bf15=robjects.vectors.FloatVector(fullhash["bf15"][pops[1]])
	f27_bf15=robjects.vectors.FloatVector(fullhash["bf15"][pops[2]])
	f37_bf15=robjects.vectors.FloatVector(fullhash["bf15"][pops[3]])
	a_bf37=[0.75,1.92,3.09,4.26]
	base_bf37=robjects.vectors.FloatVector(fullhash["bf37"][pops[0]])
	f15_bf37=robjects.vectors.FloatVector(fullhash["bf37"][pops[1]])
	f27_bf37=robjects.vectors.FloatVector(fullhash["bf37"][pops[2]])
	f37_bf37=robjects.vectors.FloatVector(fullhash["bf37"][pops[3]])
	a_f15f37=[1,2.17,3.34,4.5]
	base_f15f37=robjects.vectors.FloatVector(fullhash["f15f37"][pops[0]])
	f15_f15f37=robjects.vectors.FloatVector(fullhash["f15f37"][pops[1]])
	f27_f15f37=robjects.vectors.FloatVector(fullhash["f15f37"][pops[2]])
	f37_f15f37=robjects.vectors.FloatVector(fullhash["f15f37"][pops[3]])
	
	mean_bf15=robjects.vectors.FloatVector(newhash["bf15"]["mean"])
	mean_bf37=robjects.vectors.FloatVector(newhash["bf37"]["mean"])
	mean_f15f37=robjects.vectors.FloatVector(newhash["f15f37"]["mean"])
	stdabw_bf15=robjects.vectors.FloatVector(newhash["bf15"]["stdabw"])
	stdabw_bf37=robjects.vectors.FloatVector(newhash["bf37"]["stdabw"])
	stdabw_f15f37=robjects.vectors.FloatVector(newhash["f15f37"]["stdabw"])
	se_bf15=robjects.vectors.FloatVector(newhash["bf15"]["se"])
	se_bf37=robjects.vectors.FloatVector(newhash["bf37"]["se"])
	se_f15f37=robjects.vectors.FloatVector(newhash["f15f37"]["se"])
	#gplot.plotCI(mean_bf15,uiw=stdabw_bf15,at=robjects.vectors.FloatVector(a_bf15),cex=4,ylim=r.range(0,1),col="green",xlab=popnames)
	#gplot.plotCI(mean_bf37,uiw=stdabw_bf37,at=robjects.vectors.FloatVector(a_bf37),cex=4,col="red",add=True)
	#gplot.plotCI(mean_f15f37,uiw=stdabw_f15f37,at=robjects.vectors.FloatVector(a_f15f37),cex=4,col="blue",add=True)
	#graphics.axis(cex_axis=1.15)
	graphics.par(cex_axis = 2)
	graphics.par1(cex_lab = 2)
	
	graphics.boxplot(base_bf15,f15_bf15,f27_bf15,f37_bf15,at=robjects.vectors.FloatVector(a_bf15),col="green",boxwex=0.2,ylab="allele_frequencies [%]" )
	graphics.boxplot(base_bf37,f15_bf37,f27_bf37,f37_bf37,at=robjects.vectors.FloatVector(a_bf37),col="red",add=True,boxwex=0.2,names=popnames)
	graphics.boxplot(base_f15f37,f15_f15f37,f27_f15f37,f37_f15f37,at=robjects.vectors.FloatVector(a_f15f37),col="blue",add=True,boxwex	=0.2)
	graphics.legend(x= 2.35, y = 1,legend=candr,col=colr,cex=2,pch=15)
	grdevices.dev_off()
if len(str(options.pops).split(","))==3:
	pops=str(options.pops).split(",")
	cand,fullhash,newhash=["bf15","bf37","f15f37"],{},{}
	feat=["mean","stdabw","se"]
	for l in cand:
		fullhash[l]=dict(zip(pops,[[],[],[],[]]))
		newhash[l]=dict(zip(feat,[[],[],[]]))
	for l in open(str(options.sync).split(",")[0],"r"):
		for i in range(0,len(pops)):
			#print l.rstrip()+"\t"+str(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))])+"\t"+str(max_allele_change(l,int(pops[0]),int(pops[1])))+"\t"+str(pops[1])
			fullhash["bf15"][pops[i]].append(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))])
			#print freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))]
	for l in open(str(options.sync).split(",")[1],"r"):
		for i in range(0,len(pops)):
			#print l.rstrip()+"\t"+str(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))])+"\t"+str(max_allele_change(l,int(pops[0]),int(pops[1])))+"\t"+str(pops[1])
			fullhash["bf37"][pops[i]].append(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[2]))])
			#print fullhash
	for l in open(str(options.sync).split(",")[2],"r"):
		for i in range(0,len(pops)):
			#print l.rstrip()+"\t"+str(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[0]),int(pops[1]))])+"\t"+str(max_allele_change(l,int(pops[0]),int(pops[1])))+"\t"+str(pops[1])
			fullhash["f15f37"][pops[i]].append(freqhash(l,int(pops[i]))[max_allele_change(l,int(pops[1]),int(pops[2]))])
	#print fullhash
	for k,v in sorted(fullhash.items()):
		for k1,v1 in sorted(v.items()):
			newhash[k]["mean"].append(meanstdv(v1)[0])
	for k,v in sorted(fullhash.items()):
		for k1,v1 in sorted(v.items()):
			newhash[k]["stdabw"].append(meanstdv(v1)[1])
	for k,v in sorted(fullhash.items()):
		for k1,v1 in sorted(v.items()):
			newhash[k]["se"].append(meanstdv(v1)[2])
	#print newhash
	
	r=robjects.r
	graphics = importr('graphics')
	grdevices = importr('grDevices')
	gplot=importr('gplots')
	graphics.par = SignatureTranslatedFunction(graphics.par,init_prm_translate = {'cex_axis': 'cex.axis'})
	graphics.par1 = SignatureTranslatedFunction(graphics.par,init_prm_translate = {'cex_lab': 'cex.lab'})
	popnames = robjects.StrVector(pop)
	print popnames
	#sys.exit()
	grdevices.pdf(file=out,width=30,height=10)
	a_bf15=[0.5,1.67,2.84]
	candr=robjects.vectors.StrVector(cand)
	colr=robjects.vectors.StrVector(["green","red","blue"])
	
	base_bf15=robjects.vectors.FloatVector(fullhash["bf15"][pops[0]])
	f15_bf15=robjects.vectors.FloatVector(fullhash["bf15"][pops[1]])
	f37_bf15=robjects.vectors.FloatVector(fullhash["bf15"][pops[2]])
	a_bf37=[0.75,1.92,3.09]
	base_bf37=robjects.vectors.FloatVector(fullhash["bf37"][pops[0]])
	f15_bf37=robjects.vectors.FloatVector(fullhash["bf37"][pops[1]])
	f37_bf37=robjects.vectors.FloatVector(fullhash["bf37"][pops[2]])
	a_f15f37=[1,2.17,3.34]
	base_f15f37=robjects.vectors.FloatVector(fullhash["f15f37"][pops[0]])
	f15_f15f37=robjects.vectors.FloatVector(fullhash["f15f37"][pops[1]])
	f37_f15f37=robjects.vectors.FloatVector(fullhash["f15f37"][pops[2]])
	
	mean_bf15=robjects.vectors.FloatVector(newhash["bf15"]["mean"])
	mean_bf37=robjects.vectors.FloatVector(newhash["bf37"]["mean"])
	mean_f15f37=robjects.vectors.FloatVector(newhash["f15f37"]["mean"])
	stdabw_bf15=robjects.vectors.FloatVector(newhash["bf15"]["stdabw"])
	stdabw_bf37=robjects.vectors.FloatVector(newhash["bf37"]["stdabw"])
	stdabw_f15f37=robjects.vectors.FloatVector(newhash["f15f37"]["stdabw"])
	se_bf15=robjects.vectors.FloatVector(newhash["bf15"]["se"])
	se_bf37=robjects.vectors.FloatVector(newhash["bf37"]["se"])
	se_f15f37=robjects.vectors.FloatVector(newhash["f15f37"]["se"])
	#gplot.plotCI(mean_bf15,uiw=stdabw_bf15,at=robjects.vectors.FloatVector(a_bf15),cex=4,ylim=r.range(0,1),col="green",xlab=popnames)
	#gplot.plotCI(mean_bf37,uiw=stdabw_bf37,at=robjects.vectors.FloatVector(a_bf37),cex=4,col="red",add=True)
	#gplot.plotCI(mean_f15f37,uiw=stdabw_f15f37,at=robjects.vectors.FloatVector(a_f15f37),cex=4,col="blue",add=True)
	#graphics.axis(cex_axis=1.15)
	graphics.par(cex_axis = 2)
	graphics.par1(cex_lab = 2)
	graphics.boxplot(base_bf15,f15_bf15,f37_bf15,at=robjects.vectors.FloatVector(a_bf15),col="green",boxwex=0.2,ylab="allele_frequencies [%]" )
	graphics.boxplot(base_bf37,f15_bf37,f37_bf37,at=robjects.vectors.FloatVector(a_bf37),col="red",add=True,boxwex=0.2,names=popnames)
	graphics.boxplot(base_f15f37,f15_f15f37,f37_f15f37,at=robjects.vectors.FloatVector(a_f15f37),col="blue",add=True,boxwex	=0.2)
	graphics.legend(x= 2.35, y = 1,legend=candr,col=colr,cex=2,pch=15)
	grdevices.dev_off()