import sys
import collections
from optparse import OptionParser,OptionGroup


#Author: Martin Kapun

#########################################################   HELP   #########################################################################
#print
usage= "python %prog -m sample.mpileup -s inh.sam  -n pop1,pop2,pop3 > output.cov"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,"""
H E L P:
____________

This script calculates the average coverage of all populations in a pileup or mpileup 1) split by chromosome and 2) average across the euchromatin (2L,2R,3L,3R,X,4)

""") 
#########################################################   CODE   #########################################################################

parser.add_option("-m", "--mpileup", dest="m", help="mpileup file")
parser.add_option("-n", "--names", dest="n", help="list of names of the populations in the mpileup, need to be separated by a comma; this information will be used for the header")
parser.add_option("-s", "--samheader", dest="s", help="Header of SAM file; remember! samtools view -H")
(options, args) = parser.parse_args()



def meanstdv(x,n):
    from math import sqrt
    mean, std, se = 0, 0,0
    for a in x:
	mean = mean + a
    mean = mean / float(n)
    for a in x:
	std = std + (a - mean)**2
    std = sqrt(std / float(n-1))
    se=std/sqrt(n)
    return mean, std,se

code=["2L","2R","3L","3R","4","X"]

leng=collections.defaultdict(lambda:0)    
totallength=0
for l in open(options.s,"r"):
	if "@" in l:
		if "@SQ" in l:
			chrom=l.split()[1].split(":")[1]
			length=l.split()[2].split(":")[1]
			leng[chrom]=length
			if chrom in code:
				totallength+=int(length)
	else:
		break
###### identify number of populations ########

file=open(options.m,"r")
datalength=(len(file.readline().split())-3)/3
datarange=range(datalength)

###### sum up coverage per population per chromosme 
tl=collections.defaultdict(lambda:0)
leng2=collections.defaultdict(lambda:collections.defaultdict(lambda:0)) 
for l in open(options.m,"r"):
	a=l.split()
	if a[0] not in leng.keys():
		continue
	for i in range(len(datarange)):
		leng2[a[0]][i]+=int(a[3+(i*3)])
		if a[0] in code:
			tl[i]+=int(a[3+(i*3)])

##### print average coverage
print "chrom\t"+"\t".join(options.n.split(","))

for k,v in sorted(leng2.items()):
	ls=[]
	for pop,su in sorted(v.items()):
		ls.append(str(su/float(leng[k])))
	print k+"\t"+"\t".join(ls)
print "______________________________"

ls=[]
for pop,su in sorted(tl.items()):
	ls.append(str(su/float(totallength)))

print "totaleuchr\t"+"\t".join(ls)
