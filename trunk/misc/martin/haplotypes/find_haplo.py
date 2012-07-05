import sys
import re
import collections 
from optparse import OptionParser, OptionGroup
import time
import datetime
import copy

#version 1.0
#Author: Martin Kapun

#########################################################   HELP   #########################################################################
print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python find_haplo.py -i haplo.inp -m Base.sam -o haplo.out
2)	This script parses a SAM file line by line and extracts counts of haplotypes from the full alignment defined by the input files (-i). Haplotypes are counted, if both SNPs are on one read. Each haplotype and its corresponding count and frequency will be printed in a new line: See example:

3R	19630309	19630309	A	18	0.5
3R	19630309	19630309	G	18	0.5
3R	19347063	19347063,19347066	GA	12	0.387096774194
3R	19347063	19347063,19347066	TA	19	0.612903225806
3R	15212908	15212908,15212913	GT	23	0.560975609756
3R	15212908	15212908,15212913	AT	18	0.439024390244

col1:	chromosome
col2:	candidate SNP position
col3:	haplotype positions (if comprised by single SNP, just one position)
col4:	haplotype (if comprised by single SNP, just one nucleotide)
col5:	haplotype count
col6:	haplotype frequency
	""")
	
#########################################################   CODE   #########################################################################


parser.add_option("-i", "--inp", dest="inp", help="input file: provide output from extract_haplo.py without endings: \"_u\" and \"_d\"")
parser.add_option("-m", "--sam", dest="sam", help="sam-file")
parser.add_option("-o", "--out", dest="out", help="output file")	
parser.add_option_group(group)
(options, args) = parser.parse_args()
#HWUSI-EAS300R:3:61:1581:1623#0	163	2L	1	29	41S33M	=	103	176	TTCGTCACGTACCAGGGGAGCCCAGATGCGTGACGTGCNGCTCTCGACTGGACCACACGTAGCCTATTAAGCTG	aab`abaaa]]baaba````aa[aV^`aa`_a`a`]aTDY^aZZa`X_\aV`X]]]ZVaU^aV^\\^Z_]_XPa	XT:A:M	NM:i:1	SM:i:29	AM:i:29	XM:i:1	XO:i:0	XG:i:0	MD:Z:18G14	XA:Z:U,+20312035,74M,2;U,+14839290,74M,3;	RG:Z:s_3
#HWUSI-EAS300R:3:70:376:588#0	163	2L	1	29	7S44M	=	159	231	GTGCAGCTCTCGACTGGACCACACGTAGCCTATTAAGCTGGGATTTGGCGG	_\_a``````_`T\`V`\]_```^`XU]`]YVSY]K^\Z]\\PS\R[[R^X	XT:A:M	NM:i:1	SM:i:29	AM:i:29	XM:i:1	XO:i:0	XG:i:0	MD:Z:18G25	XA:Z:3L,+21077701,51M,1;U,+21512660,51M,1;U,-19640527,51M,1;U,+3486670,51M,2;	RG:Z:s_3
#HWUSI-EAS613-R:8:7:316:308#0	163	2L	1	29	4S44M	=	148	223	CAGCTCTCGACTGGACCACACGTAGCCTATTAAGCTGGGATTTGGCGG	`[aaaaa[aXX_a^aaaY`_`aa[XS_\P_aaaaa^S`_aZaaa`aa_	XT:A:M	NM:i:1	SM:i:29	AM:i:29	XM:i:1	XO:i:0	XG:i:0	MD:Z:18G25	XA:Z:U,+3486673,48M,1;U,+20958764,48M,1;U,+21512663,48M,1;U,-19640527,48M,1;	RG:Z:s_8


### count number of times an item occurs in a list
def countdp(l):
	h=collections.defaultdict(int)
	for x in l:
		h[x]+=1
	return h.items()

def alignment(l,seq):
	nseq=""
	for a in l: # rebuild aligned sequence based on cigar
		if a[0]=="M":
			nseq=nseq+seq[:(int(a[1]))]
			seq=seq[int(a[1]):]
		if a[0]=="S":
			seq=seq[(int(a[1])):]
		if a[0]=="D":
			nseq=nseq+(int(a[1])*"-")
		if a[0]=="I":
			seq=seq[(int(a[1])):]
	return nseq
def igvdict(inp):
	igvd=collections.defaultdict(lambda:[])
	for l in inp: # read line in SNPs fil
		chrom=l.split()[0]
		snp=l.split()[1]
		hap=l.split()[2]
		if len(hap.split(","))!=1:
			igvd[chrom].append(hap)
	return igvd
### link proximate SNPs upstream and downstream of the candidate SNP if the are within the range sys.argv[4]  
ou=open(str(options.inp)+"_u","r")
od=open(str(options.inp)+"_d","r")

igvdictu=igvdict(ou)
igvdictd=igvdict(od)

print "start parsing SAM file:"
### extract haplotypes based on SNPs in the proximity of candidate SNP and count frequency in population
count=0
outu=open(options.out+"_u","w")
outd=open(options.out+"_d","w")
sam=open(options.sam,"r") # open sam file
canddu=collections.defaultdict(lambda:[])
canddd=collections.defaultdict(lambda:[])
t3=0
t2=0
for l in sam: # loop throgh sam
	if "@" not in l.split()[0]:
		count+=1
		qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual=l.split("\t")[0:11]	
		cigarlist=[]
		i=""
		for a in cigar: # parse cigar string
			if re.match(r"\d",a):
				i+=str(a)
			else:
				tup=(a,i)
				cigarlist.append(tup)
				i=""

		nseq=alignment(cigarlist,seq)
		alrange=range(int(pos),int(pos)+len(nseq)-1)
		aldict=dict(zip(alrange,map(str,nseq)))
		haplotype=""
		
		if rname in igvdictu:
			for values in igvdictu[rname]:
				if int(max(values.split(",")))+200<int(pos):
					igvdictu[rname].remove(values)
				if len(values.split(","))==2 and int(min(values.split(","))) in aldict and int(max(values.split(","))) in aldict: #select SAM lines that are within range around candidate snps
					haplotype+=aldict[int(min(values.split(",")))]+aldict[int(max(values.split(",")))]
					canddu[rname+"%"+values.split(",")[1]+"%"+str(values)].append(haplotype)
					#print canddu
				haplotype=""
				
		if rname in igvdictd:
			for values in igvdictd[rname]:
				if int(max(values.split(",")))+200<int(pos):
					igvdictd[rname].remove(values)
				if len(values.split(","))==2 and int(min(values.split(","))) in aldict and int(max(values.split(","))) in aldict: #select SAM lines that are within range around candidate snps
					haplotype+=aldict[int(min(values.split(",")))]+aldict[int(max(values.split(",")))]
					canddd[rname+"%"+values.split(",")[0]+"%"+str(values)].append(haplotype)
				haplotype=""
		if count%100000==0:
			t1=time.clock()
			print str(count)+" lines in SAM processed; time elapsed: "+str(datetime.timedelta(seconds=t1-t2))
			t3+=t1-t2
			t2=t1
print "done! in: "+str(datetime.timedelta(seconds=t3))

for key, value in sorted(canddu.items()):
	for hap,count in countdp(value): #print output

		outu.write(key.split("%")[0]+"\t"+str(key.split("%")[1])+"\t"+key.split("%")[2]+"\t"+hap+"\t"+str(count)+"\t"+str(float(count)/len(value))+"\n")
for key, value in sorted(canddd.items()):

	for hap,count in countdp(value): #print output
		outd.write(key.split("%")[0]+"\t"+str(key.split("%")[1])+"\t"+key.split("%")[2]+"\t"+hap+"\t"+str(count)+"\t"+str(float(count)/len(value))+"\n")		
