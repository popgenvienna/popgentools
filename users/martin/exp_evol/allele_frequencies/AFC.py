
import sys
import collections
import modules.RCMH
from modules.CMH import CMHReader
from optparse import OptionParser, OptionGroup


#Author: Martin Kapun


#########################################################   HELP   #########################################################################
usage="python %prog --input candidates_BE.cmh --allpop 1,2,3,4,5,6,7,8,9,10 --start 1,2,3 --end 8,9,10 --replicate 1,5,7,9 --names Base,F15,F27,F37 > candidates_BE_rep1.afc"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,"""	

H E L P :
_________

Description:
This script caluclates the allele frequency changes between all combinations of populations for one replicate (--replicate). The user has to provide the columns in the input (--input; has to be the output of CMH test), which are corresponding to the populations used for  SNP calling (--allpop) and the replicated populations in the starting time-point (--start) and  the end (--end). Additionally you need to provide a vector of names corresponding to the timepoints compared in the replicate. Note, that SNPs are conditioned for the alleles rising from start to end
The output will be a tab delimited file consisting of columns for all time-point combinations and the corresponding allele frequency changes. Per default the script is expecting a sync file as an input. However, it also handles CMH outputs (with a P-value at the end). Then, you need to put --type CMH
""") 
#########################################################   CODE   #########################################################################

parser.add_option("-i","--input", dest="input", help="A synchronized input file")
parser.add_option("-p","--allpop", dest="all", help="allpop")
parser.add_option("-s","--start", dest="s", help="start timepoint")
parser.add_option("-e","--end", dest="e", help="end timepoint")
parser.add_option("-a","--names", dest="nam", help="names of pops compared")
parser.add_option("-r","--replicate", dest="com", help="populations to compare")
parser.add_option("--type",dest="type",help="file type: cmh output or sync file",default="sync")
parser.add_option_group(group)
(options, args) = parser.parse_args()



comp=map(int,options.com.split(","))
allpops=map(int,options.all.split(","))
names=options.nam.split(",")
start=map(int,options.s.split(","))
end=map(int,options.e.split(","))
namel=[]
for i in range(len(names)):
	for j in range(len(names)):
		if j>i:
			namel.append(names[i]+"_"+names[j])
			
print "\t".join(namel)


if options.type=="sync":
	filehandle=SyncReader(options.input)
else:
	filehandle=CMHReader(options.input)

for sync in filehandle:
	afclist=[]
	# extract the correct populations
	pops=sync.subpopulations(allpops)
	com=sync.subpopulations(comp)
	# obtain counts for the two major alleles
	(alcount,ma,mi)=modules.RCMH.Utility.getMajorAlleleCount(pops)
	# test whihc allele rises in frequency:
	s=sync.subpopulations(start)
	e=sync.subpopulations(end)
	s1,e1=0.0,0.0
	#print s
	for i in range(len(s)):
		s1+=s[i].freqh[mi]
		e1+=e[i].freqh[mi]
	if e1/len(s)>s1/len(s):
		allele=mi
	else:
		allele=ma
	#fill allelefrequency list
	for i in range(len(com)):
		for j in range(len(com)):
			if j>i:
				afclist.append(com[j].freqh[allele]-com[i].freqh[allele])
	print "\t".join(map(str,afclist))
				
