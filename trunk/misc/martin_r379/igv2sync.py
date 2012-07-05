import sys
import collections
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version 1.0

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python igv2sync.py -i SNPs.igv -s Full.sync -t igv > SNPs.sync
2)	This script prints lines from a sync (output of synchronize_pileup.pl) file if the position (first and second column) is present in the IGV file (first and third column) or in the output of Roberts snp-frequency-diff.pl (first and second column). In this case use output with ending \"_rc\")
	""")
	
	
#########################################################   CODE   #########################################################################


parser.add_option("-i", "--inp", dest="inp", help="input file")
parser.add_option("-t", "--type", dest="type", help="input type: \"igv\" or \"snp_diff\"")
parser.add_option("-s", "--syn", dest="sync", help="output from Robert's syncronization script")
parser.add_option_group(group)
(options, args) = parser.parse_args()


cand=open(str(options.inp),"r")
#3L	2580489	2580931	2	33.880213	18.148094	27.825109	2580931 	0.72222222 	0.35672515 	0.14285714 	0.16164226 	0.04193361 	0.43830475 	10 	4 	4 	7 	35 	2555775 	2586540 	FBgn0010909 	msn
#4	367905	369387	19	23.595714	2.469366	16.22554	367960 	0.6953125 	0.2269103 	0.28 	0.45275328 	0.04343656 	0.10714286 	13 	4 	6 	9 	10 	349442 	379298 	FBgn0259214 	CG42314
#X	3184863	3184914	3	13.762645	5.705856	14.313621	3184863 	0.6875 	0.3195109 	0.3498452 	0.35110504 	0.28979139 	0.02255359 	5 	4 	20 	12 	70 	3070475 	3233360 	FBgn0000479 	dnc
#2L	5349680	5349819	16	68.399101	5.445149	15.916016	5349731 	0.67080231 	0.20423288 	0.23688394 	0.06158217 	0.0749179 	0.06673461 	26 	9 	27 	23 	87 	5346237 	5365039 	FBgn0016920 	nompC
candhash={}
tup=()
if str(options.type)=="igv":
	for t in cand:
		if t.rstrip()!="":
			a=t.split()
			(chrom,pos)=(a[0],a[2])
			k=chrom+"-"+pos
			candhash[k]=1
if str(options.type)=="snp_diff":
	for t in cand:
		if t.rstrip()!="":
			a=t.split()
			(chrom,pos)=(a[0],a[1])
			k=chrom+"-"+pos
			candhash[k]=1

sync=open(str(options.sync),"r")
for l in sync:
	if l.rstrip()!="":
		k=l.split()[0]+"-"+l.split()[1]
		if k in candhash:
			print l.rstrip()