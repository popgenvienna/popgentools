######################### filtering a bam file by tags
import pysam 
import sys, re
# for fuzzy matching
import regex
import os 
import argparse
import collections
#########################################################   HELP   #########################################################################
parser = argparse.ArgumentParser(description='Filter a bam file by read barcodes giving bam files named after the base and the tag sequence. needs the oython packages pysam and regex to be installed. It allows to give a maximal number of substitutions, although that wont work for the first position of the barcode. the default value for permitted substitutions is 1.')

parser.add_argument("--in", dest="infile", help="bam file with tagged reads", required=True)
parser.add_argument("--tags", dest="tags", help="comma separated list of tags (eg: \"TGACCAAT,ACAGTGAT,GCCAATAT,CTTGTAAT\")", required=True)
parser.add_argument("--subs", dest="subs", type=int, help="maximal number of substitutions in tag sequence", default=1)

args = parser.parse_args()
infile = vars(args)['infile']
subs = vars(args)['subs']
tags= vars(args)['tags'].split(",")
# remove whitespace in tags, if there is any
tags=[ re.sub("\s","",x) for x in tags]
outfiles=[]
## index BAM file if necessary
if not os.path.exists(infile+".bai"):
	print "indexing "+options.input
	os.system("samtools index "+infile)
base_name=os.path.basename(infile)
# get rid of trailing bam
base_name=re.sub("\.bam(?=$)","",base_name)
samfile=pysam.Samfile(infile,"rb")
tag_files=collections.defaultdict()
tag_count=collections.defaultdict(int)
tag_pats=collections.defaultdict()
count=0
# open files for tag splitting
for tag in tags:
	tag_files[tag]=pysam.Samfile(base_name+"_"+tag+".bam","wb",template=samfile)
	tag_pats[tag]=regex.compile("(?:"+tag+"){s<="+str(subs)+"}")

## split BAM file
for entry in samfile.fetch(until_eof=True):
	read_tag=entry.qname.split('#')[-1]
	count+=1
	for tag in tags:
		if regex.match(tag_pats[tag],read_tag):
			tag_count[tag] += 1
			tag_files[tag].write(entry)
			break
print "reads read:\t{}".format(count)
tot_tag=0
for tag in tags:
	print "tag {}:\t{}".format(tag,tag_count[tag])
	tot_tag += tag_count[tag]
	tag_files[tag].close
print "reads not taggged:\t{}".format(count-tot_tag)
samfile.close() 
