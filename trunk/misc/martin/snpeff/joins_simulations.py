	
import sys
import collections

inputs=sys.argv[1].split("+")
simhash=collections.defaultdict(lambda:{"start":[],"end":[]})


for i in range (0, len(inputs)):
	for l in open(inputs[i],"r"):
		if l.rstrip()!="":
			a=l.split()
			simhash[str(a[0])+"_"+str(a[1])+"\t"+str(a[2])]["start"].append(a[3])
			simhash[str(a[0])+"_"+str(a[1])+"\t"+str(a[2])]["end"].append(a[4])
for k,v in sorted(simhash.items()):
	print "\t".join(k.split("_"))+"\t"+"\t".join(v["start"])+"\t"+"\t".join(v["end"])