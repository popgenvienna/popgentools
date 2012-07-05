import sys
import collections
import re
import numpy
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun
#version 1.1

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python compare_gff3.py -i test1.gff3 -j test2.gff3 -f fb -b 10
2)	numpy package needs to be installed: http://sourceforge.net/projects/numpy/files/
	""")
	
	
#########################################################   CODE   #########################################################################

parser.add_option("-i", "--input1", dest="input1", help="GFF3 file (has to be an exonerate output)")
parser.add_option("-j", "--input2", dest="input2", help="GFF3 file (either exonerate gff3 or flybase gff3: has to be defined with flag \"-f\")")
parser.add_option("-b", "--bins", dest="bins", help="number of histogram bins for summary")
parser.add_option("-g", "--genes", dest="genes", help="fasta output from extract-longest-isoform.py")
parser.add_option("-f", "--flag", dest="flag", help="flag: \"fb\" or \"e\"")
parser.add_option_group(group)
(options, args) = parser.parse_args()
def lhash(c,b):
	a=open(c,"r")
	hash1=collections.defaultdict(lambda:0)
	for l in a:
		if l.rstrip()!="" and len(l.split("\t"))==9:
			if l.split()[2]==str(b):
				reg=re.search(r"ID=.*;Name=(.*):\d+;Parent.*",l.split()[8])
				#print reg
				if reg is None:
					print "Error!!: regex search: re.search(r\"ID=(.*)\",l.split()[8].split(\";\")[0]) does not work in: \""+l.rstrip()+"\""
					sys.exit()
				else:	
					name=reg.group(1)
					length=int(l.split()[4])-int(l.split()[3])
					hash1[name]+=length
	return hash1

def featurehash(inp,feat,translist,flag):
	fhash=collections.defaultdict(lambda:[])
	list1=["chro","start","end","pos","genelength","exonn","exonl"]
	genehash,thash,nhash={},{},{}
	if flag=="fb":
		for l in open(translist,"r"):
			if ">" in l:
				genehash[l.split()[0][1:]]=dict(zip(list1,["","","","","",0,0]))
				thash[l.split()[1]]=l.split()[3]
				nhash[l.split()[1]]=l.split()[0][1:]
		for l in open(inp,"r"):
			if l.rstrip()!="" and len(l.split("\t"))==9:
				item=l.split()
				if item[2]=="gene":
					reg=re.search(r"ID=(.*)",l.split()[8].split(";")[0])
					if reg is None:
						print "Error!!: regex search: re.search(r\"ID=(.*)\",l.split()[8].split(\";\")[0]) does not work in: \""+l.rstrip()+"\""
						sys.exit()
					else:	
						name=reg.group(1)
						fhash[name].append(l)
				if item[2]==str(feat):
					for k,v in thash.items():
						if v in l.split()[8]:
							fhash[k].append(l)
		for k,v in fhash.items():
			for line in v:
				item=line.split("\t")
				if item[2]=="gene" and k in nhash:
					genehash[nhash[k]]["chro"]=item[0]
					genehash[nhash[k]]["start"]=item[3]
					genehash[nhash[k]]["end"]=item[4]
					genehash[nhash[k]]["genelength"]=str(int(item[4])-int(item[3]))
				if item[2]=="exon":
					genehash[nhash[k]]["exonn"]+=1
					genehash[nhash[k]]["exonl"]+=int(item[4])-int(item[3])
					genehash[nhash[k]]["pos"]+=str(item[3])+"."+str(item[4])+","
	if flag=="ex":
		for l in open(translist,"r"):
			if ">" in l and l.rstrip()!="":
				genehash[l.split()[0][1:]]=dict(zip(list1,["","","","","",0,0]))
		for l in open(inp,"r"):
			if l.rstrip()!="" and len(l.split("\t"))==9:
				item=l.split("\t")
				if item[2]=="gene":
					reg=re.search(r"Name=(.*)",l.split()[8].split(";")[1])
					if reg is None:
						print "Error!!: regex search: re.search(r\"ID=(.*)\",l.split()[8].split(\";\")[0]) does not work in: \""+l.rstrip()+"\""
						sys.exit()
					else:	
						name=reg.group(1)
						genehash[name]["chro"]=item[0]
						genehash[name]["start"]=item[3]
						genehash[name]["end"]=item[4]
						genehash[name]["genelength"]=str(int(item[4])-int(item[3]))
				if item[2]==str(feat):
					reg=re.search(r"Name=(.*):\d+",l.split()[8].split(";")[1])
					if reg is None:
						print "Error!!: regex search: re.search(r\"ID=(.*)\",l.split()[8].split(\";\")[0]) does not work in: \""+l.rstrip()+"\""
						sys.exit()
					else:	
						name=reg.group(1)
						genehash[name]["exonn"]+=1
						genehash[name]["exonl"]+=int(item[4])-int(item[3])
						genehash[name]["pos"]+=str(item[3])+"."+str(item[4])+","
	return genehash


def genehash(c):
	a=open(c,"r")
	hash1=collections.defaultdict(lambda:[])
	for l in a:
		if l.rstrip()!="" and len(l.split("\t"))==9:
			if l.split()[2]=="gene":
				reg=re.search(r"Name=(.*)",l.split()[8].split(";")[1])
				if reg is None:
					print "Error!!: regex search: re.search(r\"ID=(.*)\",l.split()[8].split(\";\")[0]) does not work in: \""+l.rstrip()+"\""
					sys.exit()
				else:	
					name=reg.group(1)
					length=int(l.split()[4])-int(l.split()[3])
					pos=str(l.split()[3])+","+str(l.split()[4])
					chrom=l.split()[0]
					tup=(chrom,length,pos,name)
					hash1[name].append(tup)	
	return hash1



print "input1: "+str(options.input1)
print "input2: "+str(options.input2)
bins=options.bins
namehash={}
ghash1=genehash(options.input1)
ghash2=genehash(options.input2)
ehash1=lhash(options.input1,"exon")
ehash2=lhash(options.input2,"exon")	

print
print """Genes with more than one hit of similar score:"""
print """______________________________________________________"""
print """Gene	chrom_Hit1	length_Hit1	chrom_Hit2	length_Hit2"""
print """______________________________________________________"""

genecount=0
for key, value in ghash1.items():
	genecount+=1
	if len(value)>=2:
		p=""
		for i in range(0,len(value)):
			p+="\t"+str(value[i][0])+"\t"+str(value[i][1])+"\t"+str(value[i][2])
		print key+p
print """______________________________________________________"""

print 
print """Genes in GFF"""
print """______________________________________________________"""
print str(options.input1.split("/")[-1])+"\t"+str(options.input2.split("/")[-1])+"""	percent_found"""
print """______________________________________________________"""
print str(len(ghash1))+"\t"+str(len(ghash2))+"\t"+str(float(len(ghash1))/float(len(ghash2))*100)+"%"

lengthlist=[]
for key, value in ghash1.items():
	if len(value)==1 and ghash2[key]!=[]:
		#print value[0][1], hash2[key][0][1]
		lengthlist.append(float(value[0][1])/float(ghash2[key][0][1])*100)

print 
print """proportion of length of whole gene aligned in """+str(options.input1.split("/")[-1])+""", relative to """+str(options.input2.split("/")[-1])
print """______________________________________________________"""
print """Bins	Counts"""
print """______________________________________________________"""
count,bin=numpy.histogram(lengthlist,int(bins),range=(0.0,150))
hist=zip(bin,count)
for bin,count in hist:
	print str(bin)+"\t"+str(count)

#print hash1,hash2	
#print hash1["CG15894"],hash2["CG15894"]
lengthlist=[]
print
print """genes with less than 10% overlap in """+str(options.input1.split("/")[-1])+""", relative to """+str(options.input2.split("/")[-1])
print """______________________________________________________"""
print """gene	length_gene_dataset1	length_gene_dataset2	proportion"""
print """______________________________________________________"""
for key, value in ghash1.items():
	if key in ghash2:
		if len(value)==1 and ghash2[key]!=[]:
			if float(value[0][1])/float(ghash2[key][0][1])*100 < 10.0:
				print str(key)+"\t"+ str(value[0][1])+"\t"+ str(ghash2[key][0][1])+"\t"+ str(float(value[0][1])/float(ghash2[key][0][1])*100)

lengthlist=[]
for key, value in ehash1.items():
	if key in ehash2 and ehash2[key]!=[]:
		lengthlist.append(float(value)/float(ehash2[key])*100)

#print str(cout)
print
print """proportion of length of all exons aligned in """+str(options.input1.split("/")[-1])+""", relative to """+str(options.input2.split("/")[-1])
print """______________________________________________________"""
print """Bins	Counts"""
print """______________________________________________________"""
count,bin=numpy.histogram(lengthlist,int(bins),range=(0.0,200))
hist=zip(bin,count)
for bin,count in hist:
	print str(bin)+"\t"+str(count)

print
print """exons with more than 100% overlap in """+str(options.input1.split("/")[-1])+""", relative to """+str(options.input2.split("/")[-1])
print """______________________________________________________"""
print """gene	length_gene_dataset1	length_gene_dataset2	proportion"""
print """______________________________________________________"""
for key, value in ehash1.items():
	if key in ehash2 and ehash2[key]!=[]:
		if float(value)/float(ehash2[key])*100 >= 100.0:
			print str(key)+"\t"+ str(value)+"\t"+ str(ehash2[key])+"\t"+ str(float(value)/float(ehash2[key])*100)
fhash1=featurehash(options.input1,"exon",options.genes,"ex")
fhash2=featurehash(options.input2,"exon",options.genes,options.flag)

print """gene	chrom_1	start_1	end_1	length_1	exon_length_1	No._of_exons_1	exon_coordinates_1	chrom_2	start_2	end_2	length_2	exon_length_2	No._of_exons_2	exon_coordinates_2"""
for k,v in sorted(fhash2.items()):
	print k+"\t"+str(fhash1[k]["chro"])+"\t"+str(fhash1[k]["start"])+"\t"+str(fhash1[k]["end"])+"\t"+str(fhash1[k]["genelength"])+"\t"+str(fhash1[k]["exonl"])+"\t"+str(fhash1[k]["exonn"])+"\t"+str(fhash1[k]["pos"][:-1])+"\t"+str(v["chro"])+"\t"+str(v["start"])+"\t"+str(v["end"])+"\t"+str(v["genelength"])+"\t"+str(v["exonl"])+"\t"+str(v["exonn"])+"\t"+str(v["pos"][:-1])
