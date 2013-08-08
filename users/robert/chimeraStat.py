import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import pysam


parser = OptionParser()
parser.add_option("--bam",dest="bam",help="A bam file as input")
(options, args) = parser.parse_args()

#file="/Volumes/Volume_4/analysis/pcr-non-pcr/bams/Dsim-PCRfree-rep1.filter.sort.bam"
file=options.bam
counter=collections.defaultdict(lambda:0)
samfile = pysam.Samfile( file, "rb" )
for ar in samfile:
    if not ar.is_read1:
        continue
    if ar.tid == ar.rnext:
        continue # discard those with equal id
    name1=samfile.getrname(ar.tid)
    name2=samfile.getrname(ar.rnext)
    key=tuple(sorted([name1,name2]))
    counter[key]+=1
    

items=[(k[1],k[0]) for k in counter.items()]
items=reversed(sorted(items))

for i in items:
    print(i[1],i[0])

        