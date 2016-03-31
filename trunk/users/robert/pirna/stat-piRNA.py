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
astecount=collections.defaultdict(lambda:0)

tesum=0
mirnacount=0
trnacount=0
rrnacount=0
mrnacount=0

astesum=0
asmirnacount=0
astrnacount=0
asrrnacount=0
asmrnacount=0

teld=collections.defaultdict(lambda:0)
mirnald=collections.defaultdict(lambda:0)
mrnald=collections.defaultdict(lambda:0)
trnald=collections.defaultdict(lambda:0)
rrnald=collections.defaultdict(lambda:0)

asteld=collections.defaultdict(lambda:0)
asmirnald=collections.defaultdict(lambda:0)
asmrnald=collections.defaultdict(lambda:0)
astrnald=collections.defaultdict(lambda:0)
asrrnald=collections.defaultdict(lambda:0)


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
     
     antisense=False
     if flag& 0x10:
          antisense=True
     ref=a[2]
     readlen=len(a[9])
     
     if ref.endswith("_te"):
          teseq=ref[:-3]
          teld[readlen]+=1
          tecount[teseq]+=1
          tesum+=1
          if antisense:
               astesum+=1
               astecount[teseq]+=1
               asteld[readlen]+=1
     elif ref.endswith("_miRNA"):
          mirnald[readlen]+=1
          mirnacount+=1
          if antisense:
               asmirnacount+=1
               asmirnald[readlen]+=1
     elif ref.endswith("_rRNA") or  ref.endswith("_rRNA;"):
          rrnald[readlen]+=1
          rrnacount+=1
          if antisense:
               asrrnacount+=1
               asrrnald[readlen]+=1
     elif ref.endswith("_tRNA"):
          trnald[readlen]+=1
          trnacount+=1
          if antisense:
               astrnacount+=1
               astrnald[readlen]+=1
     elif ref.endswith("_mRNA"):
          mrnald[readlen]+=1
          mrnacount+=1
          if antisense:
               asmrnacount+=1
               asmrnald[readlen]+=1
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

teprintlist=["1360","412","ACCORD","AF222049","AF418572","AF541951","BAGGINS","BLOOD","BS","BS3","BS4","Beagle","Beagle2","CIRC","DIVER2","DM06920","DM23420","DM33463","DM88","DMAURA","DMBARI1","DMBLPP","DMCOPIA","DMCR1A","DMDM11","DME010298","DME278684","DME487856","DME542581","DME9736","DMGYPF1A","DMHFL1","DMIFACA","DMIS176","DMIS297","DMLINEJA","DMMDG3","DMREPG","DMRER1DM","DMRER2DM","DMRTMGD1","DMTHB1","DMTN1731","DMTNFB","DMTOM1_LTR","DMTRDNA","DMU89994","DMW1DOC","DMZAM","DM_ROO","DOC2","DOC3","DOC4","DOC5","F","FB","FROGGER","FW2","FW3","G2","G3","G4_DM","G5A","G5_DM","G6_DM","G7","GTWIN","GYPSY10","GYPSY11","GYPSY12","GYPSY2","GYPSY3","GYPSY4","GYPSY5","GYPSY6","GYPSY7","GYPSY8","GYPSY9","HEL","HOPPER2","INE1","INVADER","INVADER2","INVADER3","INVADER4","INVADER5","INVADER6","IVK","JOCKEY2","JUAN","LOOPER1_DM","M14653","MARINER2","McCLINTOCK","OPUS","OSV","PPI251","Q","QBERT","QUASIMODO","R1-2","ROOA_LTR","ROVER","ROXELEMENT","RT1B","RT1C","S2","SPRINGER","STALKER","STALKER2","STALKER3","STALKER4","TABOR","TC1","TC1-2","TC3","TIRANT","TRANSIB1","TRANSIB2","TRANSIB3","TRANSIB4","Tinker"]

 
print "{0}\t{1}\t{2}\t{3}".format("gen","tecount",tesum,astesum)    
print "{0}\t{1}\t{2}\t{3}".format("gen","mirnacount",mirnacount,asmirnacount)
print "{0}\t{1}\t{2}\t{3}".format("gen","trnacount",trnacount,astrnacount)
print "{0}\t{1}\t{2}\t{3}".format("gen","rrnacount",rrnacount,asrrnacount)   
print "{0}\t{1}\t{2}\t{3}".format("gen","mrnacount",mrnacount,asmrnacount)   
for te in teprintlist:
     count=tecount[te]
     ascount=astecount[te]
     print "{0}\t{1}\t{2}\t{3}".format("teabundance",te,count,ascount)        
for tel in sorted(teld.keys()):
     count=teld[tel]
     ascount=asteld[tel]
     print "{0}\t{1}\t{2}\t{3}".format("te-ld",tel,count,ascount)        
for mil in sorted(mirnald.keys()):
     count=mirnald[mil]
     ascount=asmirnald[mil]
     print "{0}\t{1}\t{2}\t{3}".format("mirna-ld",mil,count,ascount)
for trl in sorted(trnald.keys()):
     count=trnald[trl]
     ascount=astrnald[trl]
     print "{0}\t{1}\t{2}\t{3}".format("trna-ld",trl,count,ascount)
for rrl in sorted(rrnald.keys()):
     count=rrnald[rrl]
     ascount=asrrnald[rrl]
     print "{0}\t{1}\t{2}\t{3}".format("rrna-ld",rrl,count,ascount)        
