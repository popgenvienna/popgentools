import sys
import re
import datetime
import time
from optparse import OptionParser, OptionGroup

#version 1.0
#Author: Martin Kapun

#########################################################   HELP   #########################################################################
print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python extract_longest_isoform_from_gff.py -g flybase.gff3 -j longest_isoforms.fasta > longest_isoform.gff
2)	"""	)																																																																																									

parser.add_option("-g", "--gff", dest="gff", help="GFF3 file (flybase gff)")
parser.add_option("-f", "--fasta", dest="fasta", help="output from extract-longest-isoform.py")
parser.add_option("-o", "--out", dest="out", help="output file")

parser.add_option_group(group)
(options, args) = parser.parse_args()

############################################# code ########################################################
out=open(options.out,"w")

transhash1,transhash2={},{}
for l in open(options.fasta,"r"):
	if ">" in l:
		transhash1[l.split()[1]]=1
		transhash2[l.split()[3]]=1
transhash3=dict(transhash1.items() + transhash2.items())		
count=0
c=0
k=0
gene=0
t3=0
t2=0
for l in open(options.gff,"r"):
		#print len(transhash1)
		if l.rstrip()!="" and len(l.split())==9:
			if l.split()[2]=="gene":
				gene=re.search(r"ID=(\w+);.*",l.split()[8]).group(1)
				if gene in transhash1:
					out.write(l.rstrip()+"\n")
					count+=1
					c=0
					if gene in transhash3:
						del transhash3[gene]
					if k in transhash2:
						del transhash2[k]
			if l.split()[2]!="exon":
				for k,v in transhash2.items():
					if k in l.split()[8] and str(l.split()[7])=="1":
						out.write(str(l.split()[0])+"\t"+str(l.split()[1])+"\t"+str(l.split()[2])+"\t"+str(l.split()[3])+"\t"+str(l.split()[4])+"\t"+str(l.split()[5])+"\t"+str(l.split()[6])+"\t2\t"+str(l.split()[8])+"\n")
						c+=1
						if k in transhash3:
							del transhash3[k]
						if gene in transhash1:
							del transhash1[gene]
					if k in l.split()[8] and str(l.split()[7])=="2": 
						out.write(str(l.split()[0])+"\t"+str(l.split()[1])+"\t"+str(l.split()[2])+"\t"+str(l.split()[3])+"\t"+str(l.split()[4])+"\t"+str(l.split()[5])+"\t"+str(l.split()[6])+"\t1\t"+str(l.split()[8])+"\n")
						c+=1
						if k in transhash3:
							del transhash3[k]
						if gene in transhash1:
							del transhash1[gene]
					if k in l.split()[8] and l.split()[7]==".":
						out.write(l.rstrip()+"\n")
						c+=1
						if k in transhash3:
							del transhash3[k]
						if gene in transhash1:
							del transhash1[gene]
					if k in l.split()[8] and l.split()[7]=="0":
						out.write(l.rstrip()+"\n")
						c+=1
						if k in transhash3:
							del transhash3[k]
						if gene in transhash1:
							del transhash1[gene]
		if count!=0 and count%2==0 and c==0:
			t1=time.clock()
			print str(count)+" genes processed in :"+str(datetime.timedelta(seconds=t1-t2))
			t3+=t1-t2
			t2=t1
		c+=1
		if len(l.split())==1:
			break
print "not found:"
for k,v in transhash3.items():
	print k