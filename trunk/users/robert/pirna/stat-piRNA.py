#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import math
import collections
import fileinput




parser = argparse.ArgumentParser(description="""           
Description
-----------
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""

Authors
-------
    Robert Kofler
""")
parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="A sam file")
parser.add_argument("--min-mq", type=int, required=False, dest="minmq", default=0, help="min mapping quality")
parser.add_argument("--max-mm", type=int, required=False, dest="maxmm", default=10, help="max mismatches")
args = parser.parse_args()
minmq=args.minmq
maxmm=args.maxmm

tecount=collections.defaultdict(lambda:0)
tesum=0
mirnacount=0
trnacount=0
rrnacount=0

teld=collections.defaultdict(lambda:0)
mirnald=collections.defaultdict(lambda:0)
trnald=collections.defaultdict(lambda:0)
rrnald=collections.defaultdict(lambda:0)

for line in args.sam:
     """
0         1         2              3    4         5    6         7      8            9                        10                  11
r1	16	M14653_te	172	70	23M	*	0	0	ATGTCGAGTTTCGTGCCGAATAA	FFFFFFFFFFFFFFFFFFBBBBB	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:23
r2	0	M14653_te	240	70	27M	*	0	0	AACAGCTGCGGAATCGCACCGAATGCT	BBBBBFFFFFBFFFFFFFFFFFFFFFF	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:27
     """
     a=line.rstrip("\n").split("\t")
     flag=int(a[1])
     if flag & 0x004 > 0:
          continue
     mq=int(a[4])
     if mq< minmq:
          continue
     
     mm=0
     tmp=a[11]
     b=tmp.split(" ")
     for bt in b:
          if bt.startswith("NM:i:"):
               mm=int(bt[5:])
     if(mm>maxmm):
          continue
     
     
     
     ref=a[2]
     readlen=len(a[9])
     
     if ref.endswith("_te"):
          teseq=ref[:-3]
          teld[readlen]+=1
          tecount[teseq]+=1
          tesum+=1
     elif ref.endswith("_miRNA"):
          mirnald[readlen]+=1
          mirnacount+=1
     elif ref.endswith("_rRNA") or  ref.endswith("_rRNA;"):
          rrnald[readlen]+=1
          rrnacount+=1
     elif ref.endswith("_tRNA"):
          trnald[readlen]+=1
          trnacount+=1
     elif  ref.endswith("_snoRNA;") or ref.endswith("_snoRNA") or ref.endswith("_snRNA;") or ref.endswith("_snRNA"):
          pass
     else:
          raise Exception("Unknown sequence end "+ ref)




#tecount=collections.defaultdict(lambda:0)
#tesum=0
#mirnacount=0
#trnacount=0
#rrnacount=0
#teld=collections.defaultdict(lambda:0)
#mirnald=collections.defaultdict(lambda:0)
#trnald=collections.defaultdict(lambda:0)
#rrnald=collections.defaultdict(lambda:0)
 
print "{0}\t{1}\t{2}".format("gen","tecount",tesum)    
print "{0}\t{1}\t{2}".format("gen","mirnacount",mirnacount)
print "{0}\t{1}\t{2}".format("gen","trnacount",trnacount)
print "{0}\t{1}\t{2}".format("gen","rrnacount",rrnacount)   
for te,count in tecount.items():
     print "{0}\t{1}\t{2}".format("teabundance",te,count)        
for tel in sorted(teld.keys()):
     count=teld[tel]
     print "{0}\t{1}\t{2}".format("te-ld",tel,count)        
for mil in sorted(mirnald.keys()):
     count=mirnald[mil]
     print "{0}\t{1}\t{2}".format("mirna-ld",mil,count)
for trl in sorted(trnald.keys()):
     count=trnald[trl]
     print "{0}\t{1}\t{2}".format("trna-ld",trl,count)
for rrl in sorted(rrnald.keys()):
     count=rrnald[rrl]
     print "{0}\t{1}\t{2}".format("rrna-ld",rrl,count)        
