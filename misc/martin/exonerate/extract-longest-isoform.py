from optparse import OptionParser
import re
import collections

usage = "usage: %prog --input proteom.fasta --type P > output.fasta "
parser = OptionParser(usage=usage)
parser.add_option("-i", "--input", dest="input", help="Fasta containing the transcriptom or proteom with all isoforms")
parser.add_option("-t", "--type", dest="type", help="Specify wether Fasta conatins DNA or Protein: type \"D\" or \"P\"")
(options, args) = parser.parse_args()



if str(options.type) == "P":
	fasta=open(str(options.input),"r")
	p=fasta.readline()
	reg=re.search(r">(FBpp\w*)\s.*loc=\w*:.*name=(.*)-P\w+?;\sparent=(FBgn\w+),(\w+)?;.*length=(\d*);.*",p)
	len1=reg.group(5)
	name=reg.group(2)
	fbgn=reg.group(3)
	if reg.group(4)!=None:
		tr=reg.group(4)
	else:
		tr="na"
	pr=reg.group(1)
	list1=[]
	namehash=collections.defaultdict(lambda:[])
	for l in fasta: 
			if ">" in l:
				namehash[name].append(tup)
				reg=re.search(r">(FBpp\w*)\s.*loc=\w*:.*name=(.*)-P\w+?;\sparent=(FBgn\w+),(\w+)?;.*length=(\d*);.*",l)
				len1=reg.group(5)
				name=reg.group(2)
				fbgn=reg.group(3)
				pr=reg.group(1)
				if reg.group(4)!=None:
					tr=reg.group(4)
				else:
					tr="na"
				list1=[]
				tup=()
			else:
				list1.append(l)
				tup=(len1,fbgn,list1,pr,tr)
				#print tup[0],tup[1],tup[3],tup[4]
	namehash[name].append(tup)
	for key, values in namehash.items():
		print ">"+key+"\t"+sorted(values, key = lambda a: a[0])[-1:][0][1]+"\t"+sorted(values, key = lambda a: a[0])[-1:][0][3]+"\t"+sorted(values, key = lambda a: a[0])[-1:][0][4]+"\t"+sorted(values, key = lambda a: a[0])[-1:][0][0]
		for length,fbgn, fasta,pr,tr in sorted(values, key = lambda a: a[0])[-1:]:
			for items in fasta:
				print items.rstrip()


if str(options.type) == "D":
	fasta=open(str(options.input),"r")
	p=fasta.readline()
	reg=re.search(r">(FB\w*)\s.*loc=.*name=(.*)-R\w+?;\sdbx.*length=(\d*);\sparent=(FBgn\w+);\srel.*",p)
	pr=reg.group(1) 
	name=reg.group(2)
	len1=reg.group(3)
	fbgn=reg.group(4)
	list1=[]
	namehash=collections.defaultdict(lambda:[])
	for l in fasta: 
			if ">" in l:
				namehash[name].append(tup)
				reg=re.search(r">(FB\w*)\s.*loc=.*name=(.*)-R\w+?;\sdbx.*length=(\d*);\sparent=(FBgn\w+);\srel.*",l)
				pr=reg.group(1)
				name=reg.group(2)
				len1=reg.group(3)
				fbgn=reg.group(4)
				list1=[]
				tup=()
			else:
				list1.append(l)
				tup=(len1,fbgn,list1,tr)

	namehash[name].append(tup)
	for key, values in namehash.items():
		print ">"+key+"\t"+sorted(values, key = lambda a: a[0])[-1:][0][1]+"\t"+sorted(values, key = lambda a: a[0])[-1:][0][3]+"\t"+sorted(values, key = lambda a: a[0])[-1:][0][0]
		for length,fbgn, fasta, pr in sorted(values, key = lambda a: a[0])[-1:]:
			for items in fasta:
				print items.rstrip()