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
        
        

def parse_line(entries,baselist,derivedlist):
        parsed=parse_sync(entries)
        
        derivedentries=[]
        baseentries=[]
        for i in derivedlist:
                key=i-1
                derivedentries.append(parsed[key])
        
        for i in baselist:
                key=i-1
                baseentries.append(parsed[key])
        return (baseentries,derivedentries)

def get_freq(pop):
        sum=float(pop[0]+pop[1])
        fmaj=float(pop[0])/sum
        fmin=float(pop[1])/sum
        return ((fmaj,fmin))

def get_het(popfreq):
        fmaj=popfreq[0]
        fmin=popfreq[1]
        assert(fmaj>=0.0 and fmaj<=1.0)
        assert(fmin>=0.0 and fmin<=1.0)
        h=1-fmaj**2-fmin**2
        return h
        
def get_fst(base,derived):
        basesum=[0,0]
        dersum=[0,0]
        
        for b in base:
                basesum[0]+=b[0]
                basesum[1]+=b[1]
        for d in derived:
                dersum[0]+=d[0]
                dersum[1]+=d[1]
        freqd=get_freq(dersum)
        freqb=get_freq(basesum)
        freqt=((freqd[0]+freqb[0])/2.0,(freqd[1]+freqb[1])/2.0)
        
        hetd=get_het(freqd)
        hetb=get_het(freqb)
        hett=get_het(freqt)
        heta=(hetb+hetd)/2.0
        
        if(hett==0):
                return 0.0
        fst=(hett-heta)/hett
        return fst


parser = OptionParser()
parser.add_option("--cmh",dest="cmh",help="A file containing the cmh results")
parser.add_option("--base",dest="base",help="A comma separated list of base populations")
parser.add_option("--derived",dest="derived",help="A comma separated list of derived populations")
(options, args) = parser.parse_args()

baselist=map(int,options.base.split(","))
derivedlist=map(int,options.derived.split(","))

for l in open(options.cmh):
        l=l.rstrip()
        a=l.split("\t")
        a.pop() #get rid of cmh-pvalue
        chr=a.pop(0)
        pos=a.pop(0)
        a.pop(0) # refchar
        baseentries,derivedentries=parse_line(a,baselist,derivedlist)
        fst=get_fst(baseentries,derivedentries)
        
        topr=[chr,pos,str(fst)]
        print "\t".join(topr)
        
        
                  
                  

  

#2L      6580    C       2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   2:0:70:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   2:0:70:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0