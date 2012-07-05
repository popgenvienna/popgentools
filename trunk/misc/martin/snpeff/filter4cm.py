import sys
import collections
candh=collections.defaultdict(lambda:"")
for l in open(sys.argv[2],"r"):
	a=l.split()
	snp=a[0]+"_"+a[1]
	candh[snp]=l

for l in open(sys.argv[1],"r"):
	a=l.split()
	if float(a[2])<1.0:
		start=int(a[1])-50000
		end=int(a[1])+50000
		for k,v in candh.items():
			chrom=k.split("_")[0]
			pos=int(k.split("_")[1])
			if chrom==a[0] and pos>=start and pos<=end:
				del candh[k]
			if chrom=="3R" and pos>=12253902 and pos <= 20565305:
				del candh[k]
for item in candh.values():
	print item.rstrip()
