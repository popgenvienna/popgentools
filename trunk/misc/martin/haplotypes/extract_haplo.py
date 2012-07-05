import sys
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

1)	usage: python extract_haplo.py -c cand.igv -s SNPs.igv -r 50 -o haplo.inp
2)	This script allows the identification of Haplotypes around candidate SNPs, it therefore parses a file with all SNPs identified (-s) and a list of candidate SNPs at the same time. Whenever a SNP is within the range (-r) to the candidate the position of both SNPs will be stored. If no SNP is sitting in the proximity, only the candidate position will be printed. The script produces two outputs: \"_u\" for Haplotypes where the candidate is sitting upstream to the next SNP and \"_d\", where the candidate SNP is sitting downstream to the next SNP. See example for a (_d) output: 	

		\"3R	999756	999756\"
		\"3R	13795765	13795730,13795765\"
		\"3R	13194499	13194495,13194499\"
		\"3R	14181778	14181776,14181778\"
		\"3R	21411182	21411161,21411182\" 
	c1:	chromsome
	c2:	candidate SNP position
	c3:	Haplotype (if adjacent SNP in range), with SNP positions separated by a \",\"
	""")

#########################################################   CODE   #########################################################################

parser.add_option("-c", "--cand", dest="cand", help="IGV file with candidate SNPs")
parser.add_option("-s", "--snp", dest="snp", help="IGV file with all SNPs")
parser.add_option("-r", "--range", dest="range", help="""allowed distance of proximate SNPs to candidate SNP""")
parser.add_option("-o", "--out", dest="out", help="output file")
parser.add_option_group(group)
(options, args) = parser.parse_args()

#HWUSI-EAS300R:3:61:1581:1623#0	163	2L	1	29	41S33M	=	103	176	TTCGTCACGTACCAGGGGAGCCCAGATGCGTGACGTGCNGCTCTCGACTGGACCACACGTAGCCTATTAAGCTG	aab`abaaa]]baaba````aa[aV^`aa`_a`a`]aTDY^aZZa`X_\aV`X]]]ZVaU^aV^\\^Z_]_XPa	XT:A:M	NM:i:1	SM:i:29	AM:i:29	XM:i:1	XO:i:0	XG:i:0	MD:Z:18G14	XA:Z:U,+20312035,74M,2;U,+14839290,74M,3;	RG:Z:s_3
#HWUSI-EAS300R:3:70:376:588#0	163	2L	1	29	7S44M	=	159	231	GTGCAGCTCTCGACTGGACCACACGTAGCCTATTAAGCTGGGATTTGGCGG	_\_a``````_`T\`V`\]_```^`XU]`]YVSY]K^\Z]\\PS\R[[R^X	XT:A:M	NM:i:1	SM:i:29	AM:i:29	XM:i:1	XO:i:0	XG:i:0	MD:Z:18G25	XA:Z:3L,+21077701,51M,1;U,+21512660,51M,1;U,-19640527,51M,1;U,+3486670,51M,2;	RG:Z:s_3
#HWUSI-EAS613-R:8:7:316:308#0	163	2L	1	29	4S44M	=	148	223	CAGCTCTCGACTGGACCACACGTAGCCTATTAAGCTGGGATTTGGCGG	`[aaaaa[aXX_a^aaaY`_`aa[XS_\P_aaaaa^S`_aZaaa`aa_	XT:A:M	NM:i:1	SM:i:29	AM:i:29	XM:i:1	XO:i:0	XG:i:0	MD:Z:18G25	XA:Z:U,+3486673,48M,1;U,+20958764,48M,1;U,+21512663,48M,1;U,-19640527,48M,1;	RG:Z:s_8
range1=int(options.range)
igv=open(options.cand,"r")

### count number of times an item occurs in a list


### link proximate SNPs upstream and downstream of the candidate SNP if the are within the range sys.argv[4]  

b=0
c=0
d=0
count=0
t2=0
t3=time.clock()
igvdictu=collections.defaultdict(lambda:[])
igvdictd=collections.defaultdict(lambda:[])
snpdictu=collections.defaultdict(lambda:{})
snpdictd=collections.defaultdict(lambda:{})
cand=open(options.cand,"r")
snps=open(options.snp,"r")
for l in snps: # read line in SNPs file
	#print l
	if "Chromosome" not in l: 
		count+=1
		b=int(l.split()[2]) #read upstream SNPs 
		l1=[]
		#print cand
		for m in cand: # read list of candidate SNPs
			if  "Chromosome" not in m and m.split()[0]==l.split()[0]:
				if int(m.split()[2])==c and b-c<=int(range1):
					igvdictd[m.split()[0]+"%"+m.split()[2]].append(str(c))
					igvdictd[m.split()[0]+"%"+m.split()[2]].append(str(b))
				if int(m.split()[2])==c and c-d<=int(range1):
					igvdictu[m.split()[0]+"%"+m.split()[2]].append(str(d))
					igvdictu[m.split()[0]+"%"+m.split()[2]].append(str(c))
				if int(m.split()[2])==c and b-c>int(range1):
					igvdictd[m.split()[0]+"%"+m.split()[2]].append(str(c))
				if int(m.split()[2])==c and c-d>int(range1):
					igvdictu[m.split()[0]+"%"+m.split()[2]].append(str(c))
				else:
					l1.append(m)

			else: 
				l1.append(m)
		cand=l1
		#print l.split()[0],l.split()[1],igvdictd
		d=c # read downstream SNP
		c=b # read central SNP
		if count%100000==0:
			t1=time.clock()
			print str(count)+" SNPs processed; time elapsed: "+str(datetime.timedelta(seconds=t1-t2))
			t3+=t1-t2
			t2=t1

print "done! in: "+str(datetime.timedelta(seconds=t3))
print "______________________________________________________________"
ou=open(str(options.out)+"_u","w")
od=open(str(options.out)+"_d","w")

for k,v in sorted(igvdictd.items()):
	if len(v)==2:
		od.write(str(k.split("%")[0])+"\t"+str(k.split("%")[1])+"\t"+str(v[0])+","+str(v[1])+"\n")
for k,v in sorted(igvdictu.items()):
	if len(v)==2:
		ou.write(str(k.split("%")[0])+"\t"+str(k.split("%")[1])+"\t"+str(v[0])+","+str(v[1])+"\n")
