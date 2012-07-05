import sys,select
import collections
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from optparse import OptionParser,OptionGroup
import copy
import time
import datetime
import math
import random


#Author: Martin Kapun
#version: 2.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python marlies_stouffer-i candidates.sync -b 100 -p 1,5,9+2,6,10+3,4,8 -s 1,2 -b 100 -t 0.1 -c 10 -m 500 -n 10 -o output
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/ 
	""") 


parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-p", "--pops", dest="pops", help="define 2 or 3 populations for which you want to calculate the allele frequency change. populations are separated by a \",\"(e.g. 1,5,9+2,6,10+3,4,8)")
parser.add_option("-s", "--se", dest="se", help="define start and end of pop comparison of CMH: e.g.: 1,2 for Base-F15, or 1,3 for base-f37")
parser.add_option("-b", "--bins", dest="bins", help="number of bins in histogram")
parser.add_option("-t", "--th", dest="th", help="MAF threshold")
parser.add_option("-c", "--count", dest="count", help="minor allele count")
parser.add_option("-m", "--maxcov", dest="maxcov", help="maximum coverage")
parser.add_option("-n", "--mincov", dest="mincov", help="minimum coverage")
parser.add_option("-o", "--out", dest="out", help="output file")


parser.add_option_group(group)
(options, args) = parser.parse_args()


#######################################################   functions   #######################################################################

def random_allele(inp,se): ### this function returns a random polymorphic allele
	h,h2,mincount=collections.defaultdict(lambda:0),collections.defaultdict(lambda:0),0
	bases=["A","T","C","G"]
	for i in se:
		h["A"]+=int(inp.split()[int(i)+2].split(":")[0])
		h["T"]+=int(inp.split()[int(i)+2].split(":")[1])
		h["C"]+=int(inp.split()[int(i)+2].split(":")[2])
		h["G"]+=int(inp.split()[int(i)+2].split(":")[3])
	cov=sum(h.values())
	for k,v in h.items():
		if v!=0.0:
			h2[k]=(v)
	if len(h2)>1:
		for k,v in h2.items():
			return k
			break
	else:
		return "na"

def minallele_count(inp,pops): ### test wether SNP is "fixed" or below minimum count threshold
	h,h2,h3,mincount,minfreq=collections.defaultdict(lambda:0),[],[],0,0
	
	for i in range(0,len(pops)):
		h["A"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[0])
		h["T"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[1])
		h["C"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[2])
		h["G"]+=int(inp.split()[int(pops[int(i)])+2].split(":")[3])
	cov=sum(h.values())
	for k,v in h.items():
		if v!=0.0:
			h2.append(v)
			h3.append(v/float(cov))
	if len(h2)>1:
		mincount=min(h2)
		minfreq=min(h3)
		return mincount,minfreq
	else:
		return "fixed"


def hist(x,bins,out,title): ##produce a histogram output: x=data, bins=Bins, out=outputfile, title=Title
	import rpy2.robjects as robjects
	from rpy2.robjects.packages import importr
	r=robjects.r
	graphics = importr('graphics')
	grdevices = importr('grDevices')
	#rang=robjects.vectors.FloatVector([-4.75,-4.74])
	values=robjects.vectors.FloatVector(x)
	grdevices.png(out+'_hist.png')
	graphics.hist(values,breaks=bins,main=title+"\n\n"+str(bins)+" bins",xlab="")
	grdevices.dev_off()


def qqplot(x,dist,out,title): ##produces a qq-plot: x= data, dist=test distribution, out, title
	import rpy2.robjects as robjects
	from rpy2.robjects.packages import importr
	r=robjects.r
	graphics = importr('graphics')
	grdevices = importr('grDevices')
	values=robjects.vectors.FloatVector(x)
	grdevices.png(out+'_qqplot_null.png')
	r.qqnorm(values,main=title)
	r.qqline(values)
	grdevices.dev_off()


def kolmogoroff(x,out): ###perform a kolmogoroff smirnov test: x=datalist
	print
	print "Kolmogoroff-Smirnov Test:"
	out1.write("Kolmogoroff-Smirnov Test:\n")
	print "_________________________"
	print
	import rpy2.robjects as robjects
	from rpy2.robjects.packages import importr
	r=robjects.r
	values=robjects.vectors.FloatVector(x)
	p= r['ks.test'](values,'pnorm',mean=0,sd=1)
	print "test value: "+str(p[0][0])
	out.write("test value: "+str(p[0][0])+"\n")
	if p[1][0]==0:
		print "p-value: < 2.2e-16"
		out.write("p-value: < 2.2e-16"+"\n")
	else: 
		print "p-value: "+str(p[1][0])
		out.write("p-value: "+str(p[1][0])+"\n")
	print "_________________________"
	out.write("_________________________"+"\n")


def meanstdv(x): ### calculate mean, median, stdev. standard error : x=datalist
	from math import sqrt
	n,mean,std,se,median = len(x),0,0,0,0
	for a in x:
		mean = mean + a
	mean = mean / float(n)
	for a in x:
		std = std + (a - mean)**2
	std = sqrt(std / float(n-1))
	se= std/sqrt(n)
	median=sorted(x)[len(x)/2]
	return mean,std,se,median
		
def freqhash(inp,i): ### returns allelefrequencies fo a given population i
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

def cov(inp,i): ##returns coverage of population i
	h=collections.defaultdict(lambda:0)
	h["A"]+=int(inp.split()[i+2].split(":")[0])
	h["T"]+=int(inp.split()[i+2].split(":")[1])
	h["C"]+=int(inp.split()[i+2].split(":")[2])
	h["G"]+=int(inp.split()[i+2].split(":")[3])
	cov=sum(h.values())
	if cov!=0:
		return cov
	else:
		return "na"

def null(inp,cp,ma,th):	 ## return list of allelefrequency differences between replicates (D_null)
	dn,maf,sen=[],[],0	
	sign=test_sign_null(inp,cp,ma)
	for k in range(0,len(cp)): ### cp: replicates within one timepoint, e.g.: [(1,2,3),(4,5,6)]; looping through timepoints
		for i in range(0,len(cp[k])): #looping through replicates and calculating allele freq differences in all possible combinations (1 and 2, 1 and 3, 2 and 3) 
			for j in range(i+1,len(cp[k])): ##
				dn.append((freqhash(inp,int(cp[k][j]))[ma]-freqhash(inp,int(cp[k][i]))[ma])*sign)
	return dn

def SNP_test(inp,populations,max_coverage,min_coverage,mincount,minfreq):
	minus,maxcov,mincov=0,0,0
	for i in range(0,len(populations)): # test if "-" is in sync and coverages for all populations below threshold
		if "-" in inp.split()[int(populations[i])+2] or "0:0:0:0:0:0" in inp.split()[int(populations[i])+2]:
			minus=1
		else:
			if cov(inp,int(populations[i]))>int(max_coverage):
				print cov(inp,int(populations[i])),int(max_coverage),int(min_coverage)
				maxcov=1
			elif cov(inp,int(populations[i]))<int(min_coverage):
				mincov=1		
	if minus!=1 and maxcov!=1 and mincov!=1:
		if minallele_count(inp,populations)!="fixed" and minallele_count(inp,populations)[0]>=int(mincount) and minallele_count(inp,populations)[1]>=float(minfreq):
			return 1
		else:
			return 0
	return 0

def test_sign(inp,ct,m):
	dr=0
	for k in range(0,len(ct)):
		dr+=freqhash(inp,int(ct[k][1]))[m]-freqhash(inp,int(ct[k][0]))[m]
	if dr<0.0:
		return -1
	else:
		return 1
def test_sign_null(inp,cp,m):
	dr=0
	for k in range(0,len(cp)): ### cp: replicates within one timepoint, e.g.: [(1,2,3),(4,5,6)]; looping through timepoints
		for i in range(0,len(cp[k])): #looping through replicates and calculating allele freq differences in all possible combinations (1 and 2, 1 and 3, 2 and 3) 
			for j in range(i+1,len(cp[k])): ##
				dr+=freqhash(inp,int(cp[k][j]))[m]-freqhash(inp,int(cp[k][i]))[m]
	if dr<0.0:
		return -1
	else:
		return 1	
#########################################################   CODE   #########################################################################


finalhash,comptime=collections.defaultdict(lambda:{}),[]
pops=str(options.pops).split("+")
start=int(str(options.se).split(",")[0])-1
end=int(str(options.se).split(",")[1])-1
o,allpops=[],[]
for line in pops:
	o.append(line.split(","))

genepops=zip(*o)
comprep=[genepops[start],genepops[end]] 																						# replicates at different timepoints
for l in comprep:	
	for i in l:
		allpops.append(i)
for l in pops:
	comptime.append((l.split(",")[start],l.split(",")[end])) 
#print allpops,comprep,comptime																	# timepoints  at different replicates
z_null,d_null_full,maf,p_null,d_null=[],[],[],[],[]
time1,time2,time3,snpcount,minus,maxcov=0,0,0,0,0,0
r=robjects.r
out1=open(str(options.out)+"_summary.txt","w")	
for l in open(str(options.inp),"r"): 																							# go line by ine through sync
	snpcount+=1
	if l.rstrip()=="":																										 	# make sure that line is not empty
		continue
	if SNP_test(l,allpops,options.maxcov,options.mincov,options.count,options.th)==1:											# test if position is a SNP
		min1=random_allele(l,allpops) 																							# define random allele 
		d_n=null(l,comprep,min1,float(options.th))																				# calculate all possible frequ differences between replicates at timepoint 1 and timepoint 2 (e.g. for 3 replicates: 6 values)
		for item in d_n:
			d_null_full.append(item)																							# list of all allele freq differences of NULL
	if snpcount%100000==0: 																										# print when 100,000 SNPs are processed
		time1=time.clock()
		print str(snpcount)+" positions processed; time elapsed: "+str(datetime.timedelta(seconds=time1-time2))
		time3+=time1-time2
		time2=time1
#d_null=random.sample(d_null_full,len(d_null_full)/2) 																			# reduce d_null values by half as demanded by Marlies
d_null=d_null_full
sd_null=meanstdv(d_null)[1] 																									# calculate the SD of D_null
se_null=meanstdv(d_null)[2] 																									# calculate the SE of D_null 
mean_null=meanstdv(d_null)[0] 																									# calculate the mean of D_null 
for delta_null in d_null:  
	z_null.append((delta_null-mean_null)/sd_null) 																				# calculate z-score for allele diffences between replicates within a timepoint
for zscore in z_null:
	p_null.append(1-r['pnorm'](zscore)[0]) 																						# calculate p-values for all z-scores
	
### lists, we have now:
#d_null: raw allele freq differences between replicates
#z_null: standard normalized raw allele freq differnces between replicates
# p_null: coresponding p- values of the z-scores
 	
#################################################### print statistics ####################################################################

hist(d_null,int(options.bins),str(options.out)+"_D_null","d-score distribution of \n"+str(options.inp).split("/")[-1])
hist(z_null,int(options.bins),str(options.out)+"_Z_null","z-score distribution of \n"+str(options.inp).split("/")[-1])
hist(p_null,int(options.bins),str(options.out)+"_p_null","p-value distribution of \n"+str(options.inp).split("/")[-1])
qqplot(z_null,"qnorm",str(options.out),"qq-plot of "+str(options.inp).split("/")[-1])
kolmogoroff(z_null,out1)
r=robjects.r
print
print "method: (delta_null-mean_null)/sd_null"
out1.write("method: (delta_null-mean_null)/sd_null"+"\n")
print "_____________________________________________"
print "z_null:mean\tz_null:median\tz_null:SD\tz_null:SE"
out1.write("_____________________________________________"+"\n")
out1.write("z_null:mean\tz_null:median\tz_null:SD\tz_null:SE"+"\n")
print str(meanstdv(z_null)[0])+"\t"+str(meanstdv(z_null)[3])+"\t"+str(meanstdv(z_null)[1])+"\t"+str(meanstdv(z_null)[2])
out1.write(str(meanstdv(z_null)[0])+"\t"+str(meanstdv(z_null)[3])+"\t"+str(meanstdv(z_null)[1])+"\t"+str(meanstdv(z_null)[2])+"\n")


#################################################### continue?? ####################################################################

print
print "\033[5;35m\033[1;35mAre you happy? Shall I continue?\n\033[0m"

print "\033[1;35mYou have 15 minutes to answer with \"yes\" or \"no\", else I'll continue anyway\033[0m"

g, o, e = select.select( [sys.stdin], [], [], 900 )

if (g):
	inps=""
	inps = sys.stdin.readline().rstrip()
	if "no" in inps:
		print "OK,bye bye :-("
		sys.exit()

snpcount=0
out2=open(str(options.out)+".stouffer","w")
r=robjects.r
print
print "\033[1;31myipee, I'll continue...\033[0m"
z_real,d_real,stouffer_pvalue,snpcount=[],[],[],0	
time1,time2,time3=0,0,0
for l in open(str(options.inp),"r"):
	z_r,stouffer,p_real=0,0,[]
	snpcount+=1
	minus,maxcov,mincov=0,0,0
	if l.rstrip()=="": 																											# make sure that line is not empty
		continue
	if SNP_test(l,allpops,options.maxcov,options.mincov,options.count,options.th)==1:
		min1=random_allele(l,allpops) 																							# define the first allele with counts as "random allele"
		fullcov,scov,cv,na=0,0.0,[],0
		for k in range(0,len(comptime)):
			#print int(comptime[k][0]),int(comptime[k][1])
			cv=[cov(l,int(comptime[k][0])),cov(l,int(comptime[k][1]))] 															# sum of mean coverage for correction factor
			if meanstdv(cv)[1]!=0.0:
				fullcov+=meanstdv(cv)[0]  
			else:
				fullcov+=float(meanstdv(cv)[0])
			cv=[]
		sign=test_sign(l,comptime,min1)
		for k in range(0,len(comptime)):
			delta_real=(freqhash(l,int(comptime[k][1]))[min1]-freqhash(l,int(comptime[k][0]))[min1])*sign 								# allele frequency difference between timepoints
			z_r=((delta_real-mean_null)/sd_null)			 																	# standard normalization of delta_real
			z_real.append(z_r)  																								# list of all z scores
			d_real.append(delta_real)																							# list of all allele frequency differences
			cv=[cov(l,int(comptime[k][0])),cov(l,int(comptime[k][1]))]
			if meanstdv(cv)[1]!=0.0:
				ni=meanstdv(cv)[0]
			else:
				ni=float(meanstdv(cv)[0]) 																						# mean coverage
			cv=[]																												
			wi=math.sqrt(ni/fullcov)																							# squareroot of ni
			stouffer+=float(wi)*z_r																								# stouffer: sum of z-scores*wi
			#stouffer+=z_r
		result=1-r['pnorm'](stouffer)[0]
		stouffer_pvalue.append(result)
		if result==0:
			out2.write(l.rstrip()+"\t"+str(1.0e-16)+"\n")
		if result==1:
			out2.write(l.rstrip()+"\t"+str(0.99999999999)+"\n")
		else:
			out2.write(l.rstrip()+"\t"+str(result)+"\n")
	if snpcount%100000==0:
		time1=time.clock()
		print str(snpcount)+" positions processed; time elapsed: "+str(datetime.timedelta(seconds=time1-time2))
		time3+=time1-time2
		time2=time1
for zsc in z_real:
	p_real.append(1-r['pnorm'](zsc)[0])		
print
print "method: (delta_real-mean_null)/sd_null"
print "_____________________________________________"
print "z_real:mean\tz_real:median\tz_real:SD\tz_real:SE"
print str(meanstdv(z_real)[0])+"\t"+str(meanstdv(z_real)[3])+"\t"+str(meanstdv(z_real)[1])+"\t"+str(meanstdv(z_real)[2])
out1.write("method: (delta_real-mean_null)/sd_null"+"\n")
out1.write("_____________________________________________"+"\n")
out1.write("z_real:mean\tz_real:median\tz_real:SD\tz_real:SE"+"\n")
out1.write(str(meanstdv(z_real)[0])+"\t"+str(meanstdv(z_real)[3])+"\t"+str(meanstdv(z_real)[1])+"\t"+str(meanstdv(z_real)[2])+"\n")
hist(z_real,int(options.bins),str(options.out)+"_Z_test","z-score distribution of \n"+str(options.inp).split("/")[-1])
hist(d_real,int(options.bins),str(options.out)+"_D_test","allele freq distribution of \n"+str(options.inp).split("/")[-1])
hist(stouffer_pvalue,int(options.bins),str(options.out)+"_stouffer_p","Stouffer's p-value distribution of \n"+str(options.inp).split("/")[-1])
hist(p_real,int(options.bins),str(options.out)+"_p_real","p-value distribution of all replicates independently of \n"+str(options.inp).split("/")[-1])
print
print "\033[1;31mI'm done, puuuuuuuuuhhhhhhhhhhh..\033[0m"
