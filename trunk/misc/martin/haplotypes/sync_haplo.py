import sys
import collections 
from optparse import OptionParser,OptionGroup
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import copy
import time
import datetime

#version 2.1
#Author: Martin Kapun

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: sync_haplo.py -a base_o,F15r4_o,F15r5_o,F23 -i cand.igv -o cand.hap
2)	required external modules: rpy2 from http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script synchronizes the haplotype counts (frequencies) of multiple populations. Pairwise comparison will be performed for all populations in the input (-a). There, all populations need to be separated by a comma without spaces. (e.g.: Base_o,F15r5_o). The script must be run for upstream and downstream haplotypes (see documentation in extract_haplo.py) separately. 
	Steps:
	a)The script prints the codes for the populations in the input
	b)As a first step, the script tests whether SNPs are missing in the haplotype datasets relative to the candidate SNP list (-i). If there are inconsistencies, the program will stop at this stage.
	c)If no SNPs are missing in the inputfiles, significant changes in haplotype frequencies will be tested by Fisher exact test for each candidate in R through the rpy2 module. The output looks as follows: 	

#Populations compared:
#pop 0: /Volumes/Untitled1/No_single_end_IndelMasked/cmh/Comparison_CMH_FST_and_FET/Base3_u
#pop 1: /Volumes/Untitled1/No_single_end_IndelMasked/cmh/Comparison_CMH_FST_and_FET/F15r4_test_u
#_______________________________________________________________
#No missing SNPs in pop 0!!!
#No missing SNPs in pop 1!!!
Chr	Pos	flag	pop_compared	dist_between_SNPs	Haplotypes	Counts_POP1	Counts_POP2	P-value
2L	19190683	1	0,1	6	TA,CA	4,31	17,20	0.00167896486019
2L	17504221	2	0,1	2	AT,AG,CG	1,1,36	0,3,20	0.188118811881
2L	17962840	2	0,1	24	CC,TC,CA	23,5,4	9,8,5	0.128712871287
2L	18565314	2	0,1	8	CA,TA,CT	32,0,2	32,1,6	0.257425742574
2L	16478423	2	0,1	2	GT,AA,GA	5,2,18	0,3,16	0.0594059405941
2L	11488527	2	0,1	5	CA,TA,CG	15,20,1	10,35,0	0.039603960396
2L	19303393	1	0,1	46	TC,AC	12,5	6,8	0.156959448968
	
lines starting with # contain (I): population codes, (II) number of missing SNPs and corresponding positions.
c1:	chromsome
c2:	candidate SNP position
c3:	flag:	0:	no adjacent SNP in the range (-r) around candidate SNP; Thus haplotype=candidate SNP (e.g: A,C)
		1:	Haplotype comprised by two SNPs, but one is invariant (e.g. AA,AC)
		2:	Haplotype comprised by two SNPs and both are variant
c4:	populations compares (e.g.: pop0 vs. pop1: 0,1)
c5:	distance between SNPs in haplotypes in basepairs
c6:	Haplotypes
c7:	counts of haplotypes in pop1
c8:	counts of haplotypes in pop2
c9:	p-value from Fisher's exact test
	""")
	
	
#########################################################   CODE   #########################################################################

parser.add_option("-a", "--hap", dest="hap", help="outputs from find_haplo.py separated by a \",\"")
parser.add_option("-i", "--inp", dest="inp", help="IGV file with candidate SNPs")
parser.add_option("-o", "--out", dest="out", help="output file")
parser.add_option_group(group)
(options, args) = parser.parse_args()
########################################## print population ID's ###########################################################################
inputs=str(options.hap).split(",")
print "#Populations compared:"
for i in range(0,len(inputs)):
	print "#pop "+str(i)+": "+str(inputs[i])
print "#_______________________________________________________________"

###################################### load candidate SNPs in a hash #######################################################################
haph=collections.defaultdict(lambda:[])
snplist={}
for l in open(options.inp,"r"):
	if "Chromosome" not in l:
		snplist[l.split()[0]+"%"+l.split()[1]]=1

############################ test whether no SNPs are missing in any populations ###########################################################
c=0
for i in range(0,len(inputs)):	
	klist=copy.deepcopy(snplist)
	for l in open(inputs[i],"r"):
		if l.split()[0]+"%"+l.split()[1] in klist.keys():
			del klist[l.split()[0]+"%"+l.split()[1]]
		else: 
			continue
	if len(klist)==0:
		print "#No missing SNPs in pop "+str(i)+"!!!"
	if len(klist)!=0:
		print "#"+str(len(klist))+" SNPs not found in pop "+str(i)+"!!!"
	for k,v in klist.items():
		print "#"+k.split("%")[0]+"\t"+k.split("%")[1]
		c+=1
	klist={}
if c!=0:
	print "!!!SNPs missing in input files!!!"	
	sys.exit()

else:
	out=open(str(options.out),"w")
	out.write( "#Populations compared:\n")
	for i in range(0,len(inputs)):
		out.write("#pop "+str(i)+": "+str(inputs[i])+"\n")
	out.write("#_______________________________________________________________\n")
	
################################## load haplotypes and distance in hash ######################################################################
	hapdist={}	
	for a,v in snplist.items():
		hp=[]
		for i in range(0,len(inputs)):		
			for l in open(inputs[i],"r"):
				if l.split()[0]+"%"+l.split()[1]==a:
					b=l.split()
					tup=(b[2],b[3])
					hp.append(tup)
		hapdist[a]=set(hp)

################################## load distance between SNPs, haplotypes and corresponding counts in hash###################################
	snphap=collections.defaultdict(lambda:[])
	hapcount={}
	for a,v in snplist.items():
		snphap=collections.defaultdict(lambda:[])
		for dist,hp in hapdist[a]:
			if "," in dist:
				distance=int(max(dist.split(",")))-int(min(dist.split(",")))
			if "," not in dist:
				distance=0
			for i in range(0,len(inputs)):
				found=0
				for l in open(inputs[i],"r"):
					if l.split()[0]+"%"+l.split()[1]==a and l.split()[3]==hp:
						snphap[str(distance)+"%"+hp].append(l.split()[4])
						found=1
				if found==0:
					snphap[str(distance)+"%"+hp].append(0)
				hapcount[a]=snphap	

######################################## calculate haplotype count and compare pairwise #######################################################	
	hap=[]		
	out.write( """Chr	Pos	flag	pop_compared	dist_between_SNPs	Haplotypes	Counts_POP1	Counts_POP2	P-value"""+"\n")	
	snpcount,time2,time3=0,0,0
	for k,v in sorted(hapcount.items()):
		snpcount+=1
		for i in range(0,len(v.values()[0])):
			for j in range(i+1,len(v.values()[0])):
				pop1,pop2=[],[]
				for a,b in v.items():
					pop1.append(b[i])
					pop2.append(b[j])
					hap.append(a.split("%")[1])
					dist=a.split("%")[0]
					
############################################################## calculate haplotype flag #######################################################	
				p1,p2={},{}
				poph1=dict(zip(pop1,hap))
				poph2=dict(zip(pop2,hap))
				if 0 in poph1:
					poph1.pop(0)
				for key,value in poph1.items():
					p1[value]=key
				if 0 in poph2:
					poph2.pop(0)
				for key,value in poph2.items():
					p2[value]=key
				pophash=dict(p1.items()+p2.items())
				snp1,snp2,haplotypelist=[],[],[]
				if len(pophash.keys()[0])==1:
					flag=0
				if len(pophash.keys()[0])>1:
					for l in set(pophash.keys()):
						snp1.append(l[0])
						snp2.append(l[1])
					if len(set(snp1))==1 or len(set(snp2))==1:
						flag=1
					else:
						flag=2 
				
				comp=pop1
				for l in pop2:
					comp.append(l)
				
################################## test whether at least two haplotypes are different in pairwis eomparison ###################################				
				fieldfill=0
				for h in range(0,len(comp)/2):
					if int(comp[h])+int(comp[h+len(comp)/2])!=0:
						fieldfill+=1

#################################################### perform FET and print ####################################################################				
				if fieldfill>1:
					r=robjects.r
					g = robjects.IntVector(comp)
					m = robjects.r['matrix'](g, nrow = 2,byrow=True)
					res = r['fisher.test'](m,conf_level = 0.95, simulate_p_value=True, B = 100)[0]
					if float(res[0])<0.01:
						res = r['fisher.test'](m,conf_level = 0.95, simulate_p_value=True, B = 10000)[0]
						if float(res[0])<0.01:
							res = r['fisher.test'](m,conf_level = 0.95, simulate_p_value=True, B = 1000000)[0]
							if float(res[0])<0.01:
								res = r['fisher.test'](m,conf_level = 0.95, simulate_p_value=True, B = 100000000)[0]
								out.write( k.split("%")[0]+"\t"+str(k.split("%")[1])+"\t"+str(flag)+"\t"+str(i)+","+str(j)+"\t"+str(dist)+"\t"+(','.join(map(str,hap)))+"\t"+(','.join(map(str,pop1[:(len(pop2))])))+"\t"+(','.join(map(str,pop2)))+"\t"+str(res[0])+"\n")
								hap=[]
							else:
								out.write( k.split("%")[0]+"\t"+str(k.split("%")[1])+"\t"+str(flag)+"\t"+str(i)+","+str(j)+"\t"+str(dist)+"\t"+(','.join(map(str,hap)))+"\t"+(','.join(map(str,pop1[:(len(pop2))])))+"\t"+(','.join(map(str,pop2)))+"\t"+str(res[0])+"\n")
								hap=[]
						else:
							out.write( k.split("%")[0]+"\t"+str(k.split("%")[1])+"\t"+str(flag)+"\t"+str(i)+","+str(j)+"\t"+str(dist)+"\t"+(','.join(map(str,hap)))+"\t"+(','.join(map(str,pop1[:(len(pop2))])))+"\t"+(','.join(map(str,pop2)))+"\t"+str(res[0])+"\n")
							hap=[]
					else:
						out.write( k.split("%")[0]+"\t"+str(k.split("%")[1])+"\t"+str(flag)+"\t"+str(i)+","+str(j)+"\t"+str(dist)+"\t"+(','.join(map(str,hap)))+"\t"+(','.join(map(str,pop1[:(len(pop2))])))+"\t"+(','.join(map(str,pop2)))+"\t"+str(res[0])+"\n")
						hap=[]
				else:
					out.write( k.split("%")[0]+"\t"+str(k.split("%")[1])+"\t"+str(flag)+"\t"+str(i)+","+str(j)+"\t"+str(dist)+"\t"+(','.join(map(str,hap)))+"\t"+(','.join(map(str,pop1[:(len(pop2))])))+"\t"+(','.join(map(str,pop2)))+"\tna\n")
					hap=[]
		if snpcount%2==0:
			time1=time.clock()
			print str(snpcount)+" candidates processed; time elapsed: "+str(datetime.timedelta(seconds=time1-time2))
			time3+=time1-time2
			time2=time1
print "Done!!! time elapsed: "+str(datetime.timedelta(seconds=time3))