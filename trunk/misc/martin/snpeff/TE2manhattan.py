import sys
import collections

def group_pval(x):
	grouplist=[]
	grouphash=collections.defaultdict(lambda:0)
	c=17
	for l in open(x,"r"):
		a=l.split()
		grouplist.append(a[4])
	for i in set(grouplist):
		grouphash[i]="1e-"+str(c)
		c+=1
	return grouphash

def allele_freq(x,up,down):
	grouphash=collections.defaultdict(lambda:0)
	for l in open(x,"r"):
		a=l.split()
		grouphash[a[-1]]=str(abs(float(a[-1])*1000)*1e-20)
	return grouphash
if sys.argv[2]=="freq":
	grouph=allele_freq(sys.argv[1],25,20)
	code={"2L":"2","2R":"3","3L":"4","3R":"5","4":"6","X":"1"}
	for l in open(sys.argv[1],"r"):
		a=l.split()
		#print abs(float(a[-1]))
		if abs(float(a[-1]))>=float(sys.argv[3]):
			print code[a[0]]+"\t"+str(a[1])+"\t"+grouph[a[-1]]+"\t"+a[4]+a[0]+str(a[1])
if sys.argv[2]=="type":
	grouph=group_pval(sys.argv[1])
	code={"2L":"2","2R":"3","3L":"4","3R":"5","4":"6","X":"1"}
	for l in open(sys.argv[1],"r"):
		a=l.split()
		#print abs(float(a[-1]))
		if abs(float(a[-1]))>=float(sys.argv[3]):
			print code[a[0]]+"\t"+str(a[1])+"\t"+grouph[a[4]]+"\t"+a[4]+a[0]+str(a[1])	      