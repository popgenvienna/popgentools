import sys
import math
from optparse import OptionParser, OptionGroup
import collections

#Author: Martin Kapun
#########################################################   HELP   #########################################################################
usage="usage: %prog  --input input.pi --window-size 1000 --length inh.sam --output input_1k.pi"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P :
_________

This script bins the results of the FST and the pi output from FST4inversions.py (--input) in non-overlapping windows of size --window-size. The total lengths of the chromsomal arms need to be provided as the header of a SAM file (--length). Because FST files contain columns with meta-information which is not necessary for binning the data, which are not present in the pi files one has to specify the input file format (--data). The script will produce two outputs. *_values contains the averaged pi or FST values, whereas *_counts contains the number of actual SNPs in the bins. 
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="pi-Data")
parser.add_option("--window-size", dest="bin", help="binsize")
parser.add_option("--length", dest="l", help="length of chromosomes")
parser.add_option("--output", dest="o", help="output files")
parser.add_option("--data", dest="d", help="data type: 'fst' or 'pi'")

parser.add_option_group(group)
(options, args) = parser.parse_args()

bin=int(options.bin)

### calculate number of bins available given the chromosome length
leng=collections.defaultdict(lambda:0)
for l in open(options.l,"r"):
	a=l.split()
	leng[a[1][3:]]=int(float(a[2][3:])/bin)+1

## identify which columns to use for either a pi or a FST file
file=open(options.input,"r")
file.readline().split()
snp=file.readline().split()
chr,pos=snp[:2]

if options.d=="fst":
	ind=(len(snp)-2)/3
	datarange=range(2,ind+2)
elif options.d=="pi":
	datarange=range(2,len(snp))

apops=collections.defaultdict(lambda:[])
countpops=collections.defaultdict(lambda:0)

out1=open(options.o+"_values.txt","w")
out2=open(options.o+"_counts.txt","w")

binhash=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:0.0))))
	
for l in file:
	a=l.split()
	chr,pos=a[:2]
	## loop through all columns defined above and store the counts and sums to later calculate the averages
	for pop in datarange:
		if a[pop]!="NAN":
			## use the divisions of positions by bins as index! Clever!! Thanx Robs
			binhash[chr][int(int(pos)/bin)][pop]["count"]+=1
			binhash[chr][int(int(pos)/bin)][pop]["value"]+=float(a[pop])
		else:
			binhash[chr][int(int(pos)/bin)][pop]["count"]+=0
			binhash[chr][int(int(pos)/bin)][pop]["value"]+=0


## now loop through the dictionary of dictionaries of dictionaries and print results
for chr,hash1 in sorted(binhash.items()):
	for i in range(0,leng[chr]):
		binlist,countlist=[],[]
		if i not in hash1:
			binlist=["0"]*len(datarange)
			countlist=["0"]*len(datarange)
		else:
			for pops,val in hash1[i].items():
				if val["count"]!=0:
					binlist.append(str(val["value"]/val["count"]))
					countlist.append(str(val["count"]))
				else:
					binlist.append("0")
					countlist.append("0")
		out1.write(chr+"\t"+str(i*bin)+"\t"+"\t".join(binlist)+"\n")
		out2.write(chr+"\t"+str(i*bin)+"\t"+"\t".join(countlist)+"\n")