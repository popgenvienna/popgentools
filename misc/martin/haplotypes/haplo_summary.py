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

#pop 0: /Volumes/Temp-1/martin/martin_haplo/Base_out_nu																														#pop 1: /Volumes/Temp-1/martin/martin_haplo/BaseR3_out_nu																														#pop 2: /Volumes/Temp-1/martin/martin_haplo/F15r4_out_nu																														#pop 3: /Volumes/Temp-1/martin/martin_haplo/F15r5_out_nu																														#pop 4: /Volumes/Temp-1/martin/martin_haplo/F15r10_out_nu																														#pop 5: /Volumes/Temp-1/martin/martin_haplo/F23r1_out_nu																														#pop 6: /Volumes/Temp-1/martin/martin_haplo/F27r5_out_nu																																																												############################################# all chromosomes ########################################################																																																													/Volumes/Temp-1/martin/martin_haplo/H_d.txt										/Volumes/Temp-1/martin/martin_haplo/H_u.txt										intersect_of_both_datastsets									comparisons	total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001		total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001		total_SNPs	f0	f1	f2	0.001		0.00001		0.0000001	0,1	281	47	44	190	3	0.015789474	2	0.010526316	2	0.010526316	281	42	63	176	3	0.017045455	2	0.011363636	2	0.011363636	142	7	11	124	2	0.016129032	2	0.016129032	2	0.0161290320,2	281	47	58	176	19	0.107954545	2	0.011363636	0	0	281	42	75	164	14	0.085365854	2	0.012195122	0	0	146	7	20	119	5	0.042016807	1	0.008403361	0	00,3	281	47	52	182	143	0.785714286	67	0.368131868	26	0.142857143	281	42	66	173	139	0.803468208	70	0.404624277	24	0.138728324	135	7	10	118	80	0.677966102	34	0.288135593	12	0.1016949150,4	281	47	44	190	11	0.057894737	4	0.021052632	1	0.005263158	281	42	57	182	13	0.071428571	3	0.016483516	1	0.005494505	148	7	12	129	4	0.031007752	2	0.015503876	1	0.0077519380,5	281	47	43	191	134	0.701570681	73	0.382198953	30	0.157068063	281	42	59	180	135	0.75	85	0.472222222	29	0.161111111	153	7	13	133	85	0.639097744	37	0.278195489	13	0.0977443610,6	281	47	54	180	63	0.35	21	0.116666667	6	0.033333333	281	42	64	175	60	0.342857143	24	0.137142857	4	0.022857143	141	7	13	121	31	0.256198347	9	0.074380165	1	0.0082644631,2	281	47	60	174	13	0.074712644	5	0.028735632	1	0.005747126	281	42	71	168	12	0.071428571	3	0.017857143	0	0	134	7	14	113	4	0.03539823	2	0.017699115	0	01,3	281	47	49	185	135	0.72972973	81	0.437837838	36	0.194594595	281	42	60	179	140	0.782122905	78	0.43575419	29	0.162011173	142	7	9	126	83	0.658730159	35	0.277777778	13	0.1031746031,4	281	47	51	183	19	0.103825137	4	0.021857923	1	0.005464481	281	42	60	179	14	0.078212291	3	0.016759777	1	0.005586592	150	7	17	126	9	0.071428571	1	0.007936508	1	0.0079365081,5	281	47	39	195	140	0.717948718	80	0.41025641	36	0.184615385	281	42	47	192	149	0.776041667	89	0.463541667	45	0.234375	149	7	5	137	91	0.664233577	44	0.321167883	13	0.0948905111,6	281	47	57	177	61	0.344632768	18	0.101694915	3	0.016949153	281	42	63	176	63	0.357954545	19	0.107954545	6	0.034090909	138	7	13	118	31	0.262711864	8	0.06779661	1	0.0084745762,3	281	47	83	151	9	0.059602649	2	0.013245033	0	0	281	42	76	163	11	0.067484663	2	0.012269939	0	0	125	7	19	99	5	0.050505051	2	0.02020202	0	02,4	281	47	70	164	19	0.115853659	2	0.012195122	0	0	281	42	71	168	14	0.083333333	3	0.017857143	0	0	134	7	21	106	5	0.047169811	1	0.009433962	0	02,5	281	47	70	164	8	0.048780488	1	0.006097561	1	0.006097561	281	42	74	165	7	0.042424242	2	0.012121212	1	0.006060606	134	7	21	106	2	0.018867925	1	0.009433962	1	0.0094339622,6	281	47	91	143	2	0.013986014	2	0.013986014	0	0	281	42	79	160	1	0.00625	1	0.00625	0	0	117	7	22	88	1	0.011363636	1	0.011363636	0	03,4	281	47	61	173	98	0.566473988	52	0.300578035	26	0.150289017	281	42	60	179	102	0.569832402	51	0.284916201	15	0.083798883	132	7	11	114	54	0.473684211	28	0.245614035	11	0.0964912283,5	281	47	60	174	8	0.045977011	1	0.005747126	1	0.005747126	281	42	62	177	9	0.050847458	2	0.011299435	1	0.005649718	141	7	16	118	5	0.042372881	1	0.008474576	1	0.0084745763,6	281	47	76	158	29	0.183544304	12	0.075949367	2	0.012658228	281	42	65	174	21	0.120689655	7	0.040229885	1	0.005747126	123	7	14	102	11	0.107843137	3	0.029411765	0	04,5	281	47	49	185	98	0.52972973	48	0.259459459	20	0.108108108	281	42	49	190	112	0.589473684	58	0.305263158	27	0.142105263	147	7	11	129	57	0.441860465	26	0.201550388	10	0.077519384,6	281	47	76	158	36	0.227848101	12	0.075949367	1	0.006329114	281	42	59	180	46	0.255555556	13	0.072222222	2	0.011111111	136	7	21	108	17	0.157407407	6	0.055555556	0	05,6	281	47	62	172	17	0.098837209	3	0.01744186	3	0.01744186	281	42	64	175	15	0.085714286	5	0.028571429	3	0.017142857	136	7	16	113	7	0.061946903	3	0.026548673	3	0.026548673

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