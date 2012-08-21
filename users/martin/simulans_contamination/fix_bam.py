######################### filtering a bam file by the ID's 

import pysam 
import sys
import os 
import collections
from optparse import OptionParser, OptionGroup

#########################################################   HELP   #########################################################################
usage="python %prog --input contaminated.bam --sim simulans_1.fastq --mel melanogaster_1.fastq"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Note that the pysam package needs to be installed (type: sudo easy_install pysam) for this. This script uses the Fwd dataset of the melanogaster Fastq file (output from Ram's script; --mel ) and the Fwd dataset of the simulans dataset (--sim) as the input and splits a BAM file accordingly in a _sim and a _mel BAM which are subsets of the original BAM. Additionally, a _missed BAM will be created which may contain reads which were not in the two fastq input files (although it is supposed to be empty). This may for expampe happen if some reads were not mapped in the gsnap approach and therefore were not in any of the fastq outputs.


""") 

parser.add_option("--input", dest="input", help="A BAM file")
parser.add_option("--mel", dest="mel", help="the melanogaster specific fwd FASTQ")
parser.add_option("--sim", dest="pops", help="the melanogaster specific fwd FASTQ")

parser.add_option_group(group)
(options, args) = parser.parse_args()

## read ID's from the fastq file
def get_ids(x):
	newhash=collections.defaultdict(lambda:"")
	fastq=open(x)
	while(True):
		id=fastq.readline()[1:-3]
		newhash[id]
		fastq.readline()
		fastq.readline()
		fastq.readline()
		if id=="":
			break
	return newhash

mel=get_ids(options.mel)
sim=get_ids(options.sim)

## index BAM file if necessary
if not os.path.exists(options.input+".bai"):
	print "indexing "+options.input
	os.system("samtools index "+options.input)

samfile=pysam.Samfile(options.input,"rb")

melout=pysam.Samfile(options.input+"_mel","wb",template=samfile)
simout=pysam.Samfile(options.input+"_sim","wb",template=samfile)
missedout=pysam.Samfile(options.input+"_missed","wb",template=samfile)
## split BAM file
for l in samfile.fetch(until_eof=True):

	if l.qname in mel:
		melout.write(l)
	elif l.qname in sim:
		simout.write(l)
	else:
		missedout.write(l)

melout.close()
simout.close()
missedout.close()
samfile.close()