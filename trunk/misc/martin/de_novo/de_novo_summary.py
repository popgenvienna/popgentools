import sys
import collections
from optparse import OptionParser,OptionGroup

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python de_novo_summary.py -i contigs_k39.fa,contigs_k41.fa -p 30,50,90
	""") 
#########################################################   CODE   #########################################################################


parser.add_option("-i", "--inp", dest="inp", help="input: one or more contig files originating from CLC_de_novo or ABySS multiple files should be separated by a \",\"")
parser.add_option("-p", "--perc", dest="perc", help="length percentiles(50=median), multiple choices should be separated by a \",\"")
parser.add_option("-m", "--min", dest="min", help="minimum contig length", default=0)

(options, args) = parser.parse_args()

def median(x,n): ### calculate percentile
    totallength=sum(x)
    length=0
    sx=sorted(x)
    while length<totallength*(float(n)/100):
	item=sx.pop()
	length+=item
    return item

def analyses(inp,nh,th):
    lh,lh2,coi=collections.defaultdict(lambda:""),collections.defaultdict(lambda:""),0
    for l in open(inp,"r"):
	if ">" in l:
	    coi=l[1:]
	else:
	    lh[coi]+=l.rstrip()
    for k,v in lh.items():
	if len(v)>=int(th):
	    lh2[k]=len(v)
    for k in nh.keys():
	nh[k]=median(lh2.values(),int(k))	
    print inp.split("/")[-1]+"\t"+str(len(lh2))+"\t"+str(sum(lh2.values()))+"\t"+str(max(lh2.values()))+"\t"+"\t".join(map(str,nh.values()))

nhash=collections.defaultdict(lambda:"")
for l in options.perc.split(","):
    nhash[l]=0
    
print "dataset\tnumber_contigs\ttotal_length\tmax_length"+"\tn"+"\tn".join(nhash.keys())

for item in options.inp.split(","):
    analyses(item,nhash,options.min)
	