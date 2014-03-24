#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from syncIO import SyncReaderMajMin


        

def parse_line(entries,baselist,derivedlist):
        parsed=entries
        
        derivedentries=[]
        baseentries=[]
        for i in derivedlist:
                key=i-1
                derivedentries.append(parsed[key])
        
        for i in baselist:
                key=i-1
                baseentries.append(parsed[key])
        return (baseentries,derivedentries)

def parse_comparestring(comp):
	a=comp.split(",")
	
	baselist=[]
	derivedlist=[]
	for e in a:
		b,d=map(int,e.split("-"))
		baselist.append(b)
		derivedlist.append(d)
	return (baselist,derivedlist)


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
        
def get_average_pwfst(base,derived):
        
        assert(len(base)==len(derived))
        fstsum=0
        pwcomparisions=0
        for i in range(0,len(base)):

                b=base[i]
                d=derived[i]
                freqb=get_freq(b)
                freqd=get_freq(d)
                freqt=((freqd[0]+freqb[0])/2.0,(freqd[1]+freqb[1])/2.0)
                        
                hetd=get_het(freqd)
                hetb=get_het(freqb)
                hett=get_het(freqt)
                heta=(hetb+hetd)/2.0
        
                if(hett!=0):
                        fst=(hett-heta)/hett
                        fstsum+=fst

                pwcomparisions+=1
        
        return fstsum/pwcomparisions


parser = OptionParser()
parser.add_option("--sync",dest="sync",help="A file containing the cmh results")
parser.add_option("--compare",dest="comp",help="A comparision string as with the cmh-test")
(options, args) = parser.parse_args()

baselist,derivedlist=parse_comparestring(options.comp)

for chr,pos,mami,s in SyncReaderMajMin(options.sync):
        baseentries,derivedentries=parse_line(s,baselist,derivedlist)
        assstat=get_average_pwfst(baseentries,derivedentries)
        
        topr=[chr,str(pos),str(assstat)]
        print "\t".join(topr)






        
        
                  
                  

  

#2L      6580    C       2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   2:0:70:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   2:0:70:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0