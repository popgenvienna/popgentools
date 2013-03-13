#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections

def parse_sync(entries):
        parsed=[]
        for e in entries:
                a=map(float,e.split(":"))
                np={'A':a[0],'T':a[1],'C':a[2],'G':a[3]}
                parsed.append(np)
        ac=0
        tc=0
        cc=0
        gc=0
        for p in parsed:
                ac+=p['A']
                tc+=p['T']
                cc+=p['C']
                gc+=p['G']
                
        tmpar=[(ac,'A'),
                (tc,'T'),
                (cc,'C'),
                (gc,'G')]
        
        tmpar=sorted(tmpar, key=lambda cs: -cs[0])
        major=tmpar[0][1]
        minor=tmpar[1][1]
        toret=[]
        for p in parsed:
                novel=(p[major],p[minor])
                toret.append(novel)
        return toret
        
        

def parse_line(entries,baseset,derivedset):
        parsed=parse_sync(entries)
        derivedentries=[]
        baseentries=[]
        for i,p in enumerate(parsed):
                key=i+1
                if key in derivedset:
                        derivedentries.append(p)
                if key in baseset:
                        baseentries.append(p)
        return (baseentries,derivedentries)
        
def get_diffstat(base,derived):
        comps=[]
        for b in base:
                bf=b[0]/(b[0]+b[1])
                for d in derived:
                        df=d[0]/(d[0]+d[1])
                        diff=df-bf
                        comps.append(diff)
        posc=0
        negc=0
        for c in comps:
                if c<0:
                        negc+=1
                if c>0:
                        posc+=1
        if(negc>0 and posc>0):
                return 0.0
        absc=[]
        for c in comps:
                absc.append(abs(c))
        
        absc=sorted(absc)
        diffstat=absc[0]
        return diffstat

parser = OptionParser()
parser.add_option("--cmh",dest="cmh",help="A file containing the cmh results")
parser.add_option("--base",dest="base",help="A comma separated list of base populations")
parser.add_option("--derived",dest="derived",help="A comma separated list of derived populations")
(options, args) = parser.parse_args()

baseset=set(map(int,options.base.split(",")))
derivedset=set(map(int,options.derived.split(",")))

for l in open(options.cmh):
        l=l.rstrip()
        a=l.split("\t")
        a.pop() #get rid of cmh-pvalue
        chr=a.pop(0)
        pos=a.pop(0)
        a.pop(0) # refchar
        baseentries,derivedentries=parse_line(a,baseset,derivedset)
        diffstat=get_diffstat(baseentries,derivedentries)
        
        topr=[chr,pos,str(diffstat)]
        print "\t".join(topr)
        
        
                  
                  

  

#2L      6580    C       2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   2:0:70:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   2:0:70:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0