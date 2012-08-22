import sys
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="\npython %prog --input input.sync --pops 1+2,3,4+5 > output.sync"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,

"""
H E L P:
____________

The purpose of this script is combine different populations in a sync file (--input). E.g. columns 1 and 2 in a sync file represent the same library sequenced two times. If you want to merge these two and keep 3,4 and 5 you need to define it the following way: --pops 1+2,3,4,5. If you for exampe want to merge 1,4,5 and keep 2 you need to put it like this: --pops 1+4+5,2. Alternatively you can also merge 1 and 2 and then 3 and 4 and keep 5. Then, it should look like this:  --pops 1+2,3+4,5
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="a sync file")
parser.add_option("--pops", dest="p", help="define populations, which should be merged (+) or just printed (,); see Help for details")


populations=[map(int,x.split("+")) for x in options.p.split(",")]

def merge_pops(line,x):
	''' merge pops in sync'''
	l1=[0,0,0,0,0,0] 
	for pop in x:
		for i in range(len(l1)):
			if line[pop+2]!="-":
				l1[i]+=int(line[pop+2].split(":")[i])
	return ":".join(map(str,l1))

for l in open(options.input,"r"):
	newpop=""
	for set in populations:
		newpop+="\t"+merge_pops(l.rstrip().split("\t"),set)
	print "\t".join(l.split("\t")[1:3])+newpop
	