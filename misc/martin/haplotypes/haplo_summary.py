import sys
import collections
from optparse import OptionParser, OptionGroup

#version 1.0
#Author: Martin Kapun

#########################################################   HELP   #########################################################################
print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python haplo_summary.py -i haplotypes_out_u -j haplotypes_out_d > haplotypes.summary
2)	This script produces summary tables for all pairwise comparisons, as well as for each contig separately. It provides the number an dfrequencies of significant SNPs with three different p-value cutoffs: 10^-3,10^-5 and 10^-7. Therefore only Haplotypes, which are informative in both SNPs are considered (see description of sync_haplo.py for a comprehensive description of the information flag). A list of pop2lations compared in the analysis will printed in front of the summary tables. See sample output of all SNPs:

#pop 0: /Volumes/Temp-1/martin/martin_haplo/Base_out_nu																														

line 1-6: description of the compared pop2lations
line 7: description of the contigs for which the following summary table is produced (all , or 2L, 2R, 3L, etc.) 
line 8:	total number of SNPs in the datasets (for -i, -j and the SNPs in both dataset that have the same flag)
line 9:	number of SNPs, where no adjacent SNP was found within the given range
line 10: number of Haplotypes where only one SNP is informative (e.g. AA,AT)
line 11: number of haplotypes, where both SNPS are informative(e.g. AA,AT,TT ). Only these Haplotypes will be used for the following analyses
line 12: header of the table comprised of the names of -i, -j and the intersect
line 13 second header showing the p-value cutoffs (0.001,0.00001  and 0.0000001) and flags (f0,f1,f2) for each dataset

line 14: table columns:
col1:	pairwise comparison
col2:	total number of candidate SNPs
col3:	haplotypes with flag 0 in -i (indicates, that no SNPs were found in the proximity of the candidate (see description in sync_haplo.py))
col4:	haplotypes with flag 1 in -i (indicates, that the second SNP in the haplotype is invariable (see description in sync_haplo.py))
col5:	haplotypes with flag 2 in -i (indicates, that the second SNP in the haplotype is invariable (see description in sync_haplo.py))
col6:	counts of haplotypes being significant with p-value <0.001 in -i
col7: 	frequency of haplotypes being significant with p-value <0.001 in -i
col8:	counts of haplotypes being significant with p-value <0.00001 in -i
col9: 	frequency of haplotypes being significant with p-value <0.00001 in -i
col10:	counts of haplotypes being significant with p-value <0.0000001 in -i
col11: 	frequency of haplotypes being significant with p-value <0.0000001 in -i
col12:	haplotypes with flag 0 in -j (indicates, that no SNPs were found in the proximity of the candidate (see description in sync_haplo.py))
col13:	haplotypes with flag 1 in -j (indicates, that the second SNP in the haplotype is invariable (see description in sync_haplo.py))
col14:	haplotypes with flag 2 in -j (indicates, that the second SNP in the haplotype is invariable (see description in sync_haplo.py))
col15:	counts of haplotypes being significant with p-value <0.001 in -j
col16: 	frequency of haplotypes being significant with p-value <0.001 in -j
col17:	counts of haplotypes being significant with p-value <0.00001 in -j
col18: 	frequency of haplotypes being significant with p-value <0.00001 in -j
col19:	counts of haplotypes being significant with p-value <0.0000001 in -j
col20: 	frequency of haplotypes being significant with p-value <0.0000001 in -j
col21:	haplotypes with flag 0 in intersect of -i and -j (same flag in both datasets) (indicates, that no SNPs were found in the proximity of the candidate (see description in sync_haplo.py))
col22:	haplotypes with flag 1 in intersect of -i and -j (same flag in both datasets) (indicates, that the second SNP in the haplotype is invariable (see description in sync_haplo.py))
col23:	haplotypes with flag 2 in intersect of -i and -j (same flag in both datasets) (indicates, that the second SNP in the haplotype is invariable (see description in sync_haplo.py))
col24:	counts of haplotypes being significant with p-value <0.001 in intersect of -i and -j
col25: 	frequency of haplotypes being significant with p-value <0.001 in intersect of -i and -j
col26:	counts of haplotypes being significant with p-value <0.00001 in intersect of -i and -j
col27: 	frequency of haplotypes being significant with p-value <0.00001 in intersect of -i and -j
col28:	counts of haplotypes being significant with p-value <0.0000001 in intersect of -i and -j
col29: 	frequency of haplotypes being significant with p-value <0.0000001 in intersect of -i and -j
	""")
	
#########################################################   CODE   #########################################################################


#2LHet	187304	1	0,1	16	TC,CC,AC	15,54,0	20,64,0	0.821782178218
#2LHet	187304	1	0,2	16	TC,CC,AC	15,54,0	0,27,0	0.029702970297
#2LHet	187304	1	0,3	16	TC,CC,AC	15,54,0	1,75,1	6.83699993163e-05
#2LHet	187304	1	0,4	16	TC,CC,AC	15,54,0	7,64,0	0.0693069306931
#2LHet	187304	1	0,5	16	TC,CC,AC	15,54,0	0,71,0	9.8699999013e-06
#2LHet	187304	1	0,6	16	TC,CC,AC	15,54,0	7,70,0	0.039603960396
#2LHet	187304	1	1,2	16	TC,CC,AC	20,64,0	0,27,0	0.00305647996944
#2LHet	187304	1	1,3	16	TC,CC,AC	20,64,0	1,75,1	9.1899999081e-06
#2LHet	187304	1	1,4	16	TC,CC,AC	20,64,0	7,64,0	0.029702970297
#2LHet	187304	1	1,5	16	TC,CC,AC	20,64,0	0,71,0	2.2799999772e-06
#2LHet	187304	1	1,6	16	TC,CC,AC	20,64,0	7,70,0	0.019801980198
#2LHet	187304	1	2,3	16	TC,CC,AC	0,27,0	1,75,1	1.0
#2LHet	187304	1	2,4	16	TC,CC,AC	0,27,0	7,64,0	0.178217821782
#2LHet	187304	1	2,5	16	TC,CC,AC	0,27,0	0,71,0	na
#2LHet	187304	1	2,6	16	TC,CC,AC	0,27,0	7,70,0	0.158415841584
#2LHet	187304	1	3,4	16	TC,CC,AC	1,75,1	7,64,0	0.039603960396
#2LHet	187304	1	3,5	16	TC,CC,AC	1,75,1	0,71,0	1.0
#2LHet	187304	1	3,6	16	TC,CC,AC	1,75,1	7,70,0	0.0792079207921
#2LHet	187304	1	4,5	16	TC,CC,AC	7,64,0	0,71,0	0.039603960396
#2LHet	187304	1	4,6	16	TC,CC,AC	7,64,0	7,70,0	1.0
#2LHet	187304	1	5,6	16	TC,CC,AC	0,71,0	7,70,0	0.0594059405941

parser.add_option("-i", "--inp1", dest="inp1", help="input file1: provide output from sync_haplo.py. e.g. with ending \"_u\"")
parser.add_option("-j", "--inp2", dest="inp2", help="input file2: provide output from sync_haplo.py. e.g. with ending \"_d\"")
parser.add_option_group(group)
(options, args) = parser.parse_args()

#########################################################   print pop2lation IDs and load first dataset in hash   #########################
data1={}
popcount=0
for l in open(str(options.inp1),"r"):
	if "#pop" in l:
		print l.rstrip()
		popcount+=1
	if "#" not in l and "Chr" not in l:
		ids=l.split()[0]+"_"+l.split()[1]+"_"+l.split()[3]
		data1[ids]=l.rstrip()
allpops=[]
#########################################################   count number of pairwise comparisons  ##########################################
for i in range(0,popcount):
	for j in range(i+1,popcount):
		allpops.append(str(i)+","+str(j))

#########################################################   prepare hashes and fill them with 0   ##########################################
count=0
flago=["0","1","2"]
flaghash2a,flaghash1a,flaghashba,aflaghash2a,aflaghash1a,aflaghashba=collections.defaultdict(lambda:{}),collections.defaultdict(lambda:{}),collections.defaultdict(lambda:{}),collections.defaultdict(lambda:{}),collections.defaultdict(lambda:{}),collections.defaultdict(lambda:{})
pop,pop23,pop25,pop27,pop13,pop15,pop17,bothpop3,bothpop5,bothpop7,fullpop={},{},{},{},{},{},{},{},{},{},{}
aflags,aflaghash2,aflaghash1,aflaghashb,apop,apop23,apop25,apop27,apop13,apop15,apop17,abothpop3,abothpop5,abothpop7={},{},{},{},{},{},{},{},{},{},{},{},{},{}
for a in allpops:
	pop23[a],pop25[a],pop27[a],pop13[a],pop15[a],pop17[a],bothpop3[a],bothpop5[a],bothpop7[a],flaghash2a[a],flaghash1a[a],flaghashba[a],aflaghash2a[a],aflaghash1a[a],aflaghashba[a]=0,0,0,0,0,0,0,0,0,dict(zip(flago,[0,0,0])),dict(zip(flago,[0,0,0])),dict(zip(flago,[0,0,0])),dict(zip(flago,[0,0,0])),dict(zip(flago,[0,0,0])),dict(zip(flago,[0,0,0]))
	apop23[a],apop25[a],apop27[a],apop13[a],apop15[a],apop17[a],abothpop3[a],abothpop5[a],abothpop7[a]=0,0,0,0,0,0,0,0,0
	
#########################################################   identify first contig in table  ################################################
chro=""
for l in open(str(options.inp2),"r"):
	if "#" not in l and "Chr" not in l:
		chro=l.split()[0]
		break
#########################################################  fill hashes line by line ########################################################
for l in open(str(options.inp2),"r"):
	if "#" not in l and "Chr" not in l:
		ids=l.split()[0]+"_"+l.split()[1]+"_"+l.split()[3]
		if chro==l.split()[0]:
			chro,pos,flag,popc,a,a,a,a,pval=l.split()
			### how many entries for each flag?
			flaghash2a[popc][str(flag)]+=1
			aflaghash2a[popc][str(flag)]+=1
			flaghash1a[data1[ids].split()[3]][data1[ids].split()[2]]+=1
			aflaghash1a[data1[ids].split()[3]][data1[ids].split()[2]]+=1
			if data1[ids].split()[2]==flag:
				flaghashba[popc][str(flag)]+=1
				aflaghashba[popc][str(flag)]+=1
			### how many entries with p-values < 0.00000001 with flag==2? 
			if "na" not in l.split()[8] and l.split()[2]=="2" and float(pval)<0.001:
				pop23[popc]+=1
				apop23[popc]+=1
			if "na" not in l.split()[8] and l.split()[2]=="2" and float(pval)<0.00001:
				pop25[popc]+=1
				apop25[popc]+=1
			if "na" not in l.split()[8] and l.split()[2]=="2" and float(pval)<0.0000001:
				pop27[popc]+=1
				apop27[popc]+=1	
				
				#print data1[ids].split()[8]
			if  "na" not in data1[ids].split()[8] and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.001:
				pop13[popc]+=1
				apop13[popc]+=1
				#print pop13
			if  "na" not in data1[ids].split()[8] and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.00001:
				pop15[popc]+=1
				apop15[popc]+=1
			if  "na" not in data1[ids].split()[8] and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.0000001:
				pop17[popc]+=1
				apop17[popc]+=1
			if  "na" not in data1[ids].split()[8] and "na" not in l.split()[8] and l.split()[2]=="2" and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.001 and float(pval)<0.001:
				bothpop3[popc]+=1
				abothpop3[popc]+=1
			if  "na" not in data1[ids].split()[8] and "na" not in l.split()[8] and l.split()[2]=="2" and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.00001 and float(pval)<0.00001:
				bothpop5[popc]+=1
				abothpop5[popc]+=1
			if  "na" not in data1[ids].split()[8] and "na" not in l.split()[8] and l.split()[2]=="2" and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.0000001 and float(pval)<0.0000001:
				bothpop7[popc]+=1
				abothpop7[popc]+=1
		else: 
			fullpop[chro]=[flaghash2a,flaghash1a,flaghashba,pop23,pop25,pop27,pop13,pop15,pop17,bothpop3,bothpop5,bothpop7]
			flaghash2a,flaghash1a,flaghashba=collections.defaultdict(lambda:{}),collections.defaultdict(lambda:{}),collections.defaultdict(lambda:{})
			pop23,pop25,pop27,pop13,pop15,pop17,bothpop3,bothpop5,bothpop7={},{},{},{},{},{},{},{},{}
			for a in allpops:
				pop23[a],pop25[a],pop27[a],pop13[a],pop15[a],pop17[a],bothpop3[a],bothpop5[a],bothpop7[a],flaghash2a[a],flaghash1a[a],flaghashba[a]=0,0,0,0,0,0,0,0,0,dict(zip(flago,[0,0,0])),dict(zip(flago,[0,0,0])),dict(zip(flago,[0,0,0]))
			chro=l.split()[0]
			chro,pos,flag,popc,a,a,a,a,pval=l.split()
			### how many entries for each flag?
			flaghash2a[popc][str(flag)]+=1
			aflaghash2a[popc][str(flag)]+=1
			flaghash1a[data1[ids].split()[3]][data1[ids].split()[2]]+=1
			aflaghash1a[data1[ids].split()[3]][data1[ids].split()[2]]+=1
			if data1[ids].split()[3]==flag:
				flaghashba[popc][str(flag)]+=1
				aflaghashba[popc][str(flag)]+=1
			### how many entries with p-values < 0.00000001 with flag==2? 
			if "na" not in l.split()[8] and l.split()[2]=="2" and float(pval)<0.001:
				pop23[popc]+=1
				apop23[popc]+=1
			if "na" not in l.split()[8] and l.split()[2]=="2" and float(pval)<0.00001:
				pop25[popc]+=1
				apop25[popc]+=1
			if "na" not in l.split()[8] and l.split()[2]=="2" and float(pval)<0.0000001:
				pop27[popc]+=1
				apop27[popc]+=1	
				
				#print data1[ids].split()[8]
			if  "na" not in data1[ids].split()[8] and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.001:
				pop13[popc]+=1
				apop13[popc]+=1
				#print pop13
			if  "na" not in data1[ids].split()[8] and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.00001:
				pop15[popc]+=1
				apop15[popc]+=1
			if  "na" not in data1[ids].split()[8] and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.0000001:
				pop17[popc]+=1
				apop17[popc]+=1
			if  "na" not in data1[ids].split()[8] and "na" not in l.split()[8] and l.split()[2]=="2" and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.001 and float(pval)<0.001:
				bothpop3[popc]+=1
				abothpop3[popc]+=1
			if  "na" not in data1[ids].split()[8] and "na" not in l.split()[8] and l.split()[2]=="2" and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.00001 and float(pval)<0.00001:
				bothpop5[popc]+=1
				abothpop5[popc]+=1
			if  "na" not in data1[ids].split()[8] and "na" not in l.split()[8] and l.split()[2]=="2" and data1[ids].split()[2]=="2" and float(data1[ids].split()[8])<0.0000001 and float(pval)<0.0000001:
				bothpop7[popc]+=1
				abothpop7[popc]+=1

fullpop[chro]=[flaghash2a,flaghash1a,flaghashba,pop23,pop25,pop27,pop13,pop15,pop17,bothpop3,bothpop5,bothpop7]
fullpop["all"]=[aflaghash2a,aflaghash1a,aflaghashba,apop23,apop25,apop27,apop13,apop15,apop17,abothpop3,abothpop5,abothpop7]
#print flaghash2a
#print fullpop
flags,pop=[],[]
flags=[fullpop["all"][0],fullpop["all"][1],fullpop["all"][2]]
pop=[fullpop["all"][3],fullpop["all"][4],fullpop["all"][5],fullpop["all"][6],fullpop["all"][7],fullpop["all"][8],fullpop["all"][9],fullpop["all"][10],fullpop["all"][11]]
print
print "############################################# all chromosomes ########################################################"
print 
print "\t"+str(options.inp2)+"\t\t\t\t\t\t\t\t\t\t"+str(options.inp1)+"\t\t\t\t\t\t\t\t\t\tintersect_of_both_datastsets"
print """comparisons	total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001		total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001		total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001"""

#print flags
for a in allpops:
	#print a
	print str(a)+"\t"+str(flags[0][a]["0"]+flags[0][a]["1"]+flags[0][a]["2"])+"\t"+str(flags[0][a]["0"])+"\t"+str(flags[0][a]["1"])+"\t"+str(flags[0][a]["2"])+"\t"+str(pop[0][a])+"\t"+str(float(pop[0][a])/flags[0][a]["2"])+"\t"+str(pop[1][a])+"\t"+str(float(pop[1][a])/flags[0][a]["2"])+"\t"+str(pop[2][a])+"\t"+str(float(pop[2][a])/flags[0][a]["2"])+"\t"+str(flags[1][a]["0"]+flags[1][a]["1"]+flags[1][a]["2"])+"\t"+str(flags[1][a]["0"])+"\t"+str(flags[1][a]["1"])+"\t"+str(flags[1][a]["2"])+"\t"+str(pop[3][a])+"\t"+str(float(pop[3][a])/flags[1][a]["2"])+"\t"+str(pop[4][a])+"\t"+str(float(pop[4][a])/flags[1][a]["2"])+"\t"+str(pop[5][a])+"\t"+str(float(pop[5][a])/flags[1][a]["2"])+"\t"+str(flags[2][a]["0"]+flags[2][a]["1"]+flags[2][a]["2"])+"\t"+str(flags[2][a]["0"])+"\t"+str(flags[2][a]["1"])+"\t"+str(flags[2][a]["2"])+"\t"+str(pop[6][a])+"\t"+str(float(pop[6][a])/flags[2][a]["2"])+"\t"+str(pop[7][a])+"\t"+str(float(pop[7][a])/flags[2][a]["2"])+"\t"+str(pop[8][a])+"\t"+str(float(pop[8][a])/flags[2][a]["2"])
for k, v in sorted(fullpop.items()):
	if k!="all":
		flags=[v[0],v[1],v[2]]
		pop=[v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11]]
		print
		print "############################################# chromosome: "+str(k)+" ########################################################"
		print 
		print "\t"+str(options.inp2)+"\t\t\t\t\t\t\t\t\t\t"+str(options.inp1)+"\t\t\t\t\t\t\t\t\t\tintersect_of_both_datastsets"
		print """comparisons	total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001		total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001		total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001"""

		#print flags
		for a in allpops:
			if flags[0][a]["2"]!=0 and flags[1][a]["2"]!=0 and flags[0][a]["2"]!=0:
				#print a
				print str(a)+"\t"+str(flags[0][a]["0"]+flags[0][a]["1"]+flags[0][a]["2"])+"\t"+str(flags[0][a]["0"])+"\t"+str(flags[0][a]["1"])+"\t"+str(flags[0][a]["2"])+"\t"+str(pop[0][a])+"\t"+str(float(pop[0][a])/flags[0][a]["2"])+"\t"+str(pop[1][a])+"\t"+str(float(pop[1][a])/flags[0][a]["2"])+"\t"+str(pop[2][a])+"\t"+str(float(pop[2][a])/flags[0][a]["2"])+"\t"+str(flags[1][a]["0"]+flags[1][a]["1"]+flags[1][a]["2"])+"\t"+str(flags[1][a]["0"])+"\t"+str(flags[1][a]["1"])+"\t"+str(flags[1][a]["2"])+"\t"+str(pop[3][a])+"\t"+str(float(pop[3][a])/flags[1][a]["2"])+"\t"+str(pop[4][a])+"\t"+str(float(pop[4][a])/flags[1][a]["2"])+"\t"+str(pop[5][a])+"\t"+str(float(pop[5][a])/flags[1][a]["2"])+"\t"+str(flags[2][a]["0"]+flags[2][a]["1"]+flags[2][a]["2"])+"\t"+str(flags[2][a]["0"])+"\t"+str(flags[2][a]["1"])+"\t"+str(flags[2][a]["2"])+"\t"+str(pop[6][a])+"\t"+str(float(pop[6][a])/flags[2][a]["2"])+"\t"+str(pop[7][a])+"\t"+str(float(pop[7][a])/flags[2][a]["2"])+"\t"+str(pop[8][a])+"\t"+str(float(pop[8][a])/flags[2][a]["2"])
			else:
				print str(a)+"\t"+str(flags[0][a]["0"]+flags[0][a]["1"]+flags[0][a]["2"])+"\t"+str(flags[0][a]["0"])+"\t"+str(flags[0][a]["1"])+"\t"+str(flags[0][a]["2"])+"\t"+str(pop[0][a])+"\tna\t"+str(pop[1][a])+"\tna\t"+str(pop[2][a])+"\tna\t"+str(flags[1][a]["0"]+flags[1][a]["1"]+flags[1][a]["2"])+"\t"+str(flags[1][a]["0"])+"\t"+str(flags[1][a]["1"])+"\t"+str(flags[1][a]["2"])+"\t"+str(pop[3][a])+"\tna\t"+str(pop[4][a])+"\tna\t"+str(pop[5][a])+"\tna\t"+str(flags[2][a]["0"]+flags[2][a]["1"]+flags[2][a]["2"])+"\t"+str(flags[2][a]["0"])+"\t"+str(flags[2][a]["1"])+"\t"+str(flags[2][a]["2"])+"\t"+str(pop[6][a])+"\tna\t"+str(pop[7][a])+"\tna\t"+str(pop[8][a])+"\tna"