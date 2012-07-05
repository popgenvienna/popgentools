import sys
import collections
from modules.CMH import PopIO
from optparse import OptionParser, OptionGroup
from optparse import OptionParser,OptionGroup
import copy

#Authors:
# Martin Kapun
# Robert Kofler

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python filter_cmh.py -y candidates.cmh -p 1,2,3 -t 2
2)	Filters a cmh file and removes snps which only occur in a single population out of the populations provided by the user
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-y", "--sync", dest="sync", help="cmh output file")
parser.add_option("-p", "--pops", dest="pops", help="define populations, separate index of the replicates with ',' e.g. 3,4,5")
parser.add_option("-t", "--th", dest="th", help="min count threshold")
parser.add_option("-o", "--out", dest="out", help="output file")
parser.add_option("--test",action="store_true", dest="test",help="run the doctest")

##
##
##

parser.add_option_group(group)
(options, args) = parser.parse_args()



def count_snps(a,pops, mincount):
	# filter empty lines
	"""
	"""
	count=0
	for i in pops:
		activepop=a[int(i)-1]
		issnp=activepop.issnp(mincount)
		count+=issnp
	return count


if(options.test):
	import doctest
	doctest.testmod(verbose=1)
	sys.exit()


pops=str(options.pops).split(",")
out_d=open(str(options.out)+"_discarded","w")
out_f=open(str(options.out)+"_filtered","w")
snpc=0
for l in open(str(options.sync),"r"):
	if l.rstrip()=="":
		continue
	snpc+=1
	p=PopIO.parse_cmhline(l)

	count=count_snps(p.populations,pops, int(options.th))	
	if count>1:
		out_f.write(l)
	else: 
		out_d.write(l)
	
	if snpc%10000==0:
		print str(snpc)+" SNPs processed"
