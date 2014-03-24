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
        
def get_associastat(base,derived):
        assert(len(base)==len(derived))
        
        associastat=0
        for i,b in enumerate(base):
                d=derived[i]
                bf=b[0]/(b[0]+b[1])
                df=d[0]/(d[0]+d[1])
                diff=df-bf
                associastat+=diff
        associastat=abs(associastat)
        return associastat

def parse_comparestring(comp):
	a=comp.split(",")
	
	baselist=[]
	derivedlist=[]
	for e in a:
		b,d=map(int,e.split("-"))
		baselist.append(b)
		derivedlist.append(d)
	return (baselist,derivedlist)
		

parser = OptionParser()
parser.add_option("--sync",dest="sync",help="A file containing the cmh results")
parser.add_option("--compare",dest="comp",help="A comparision string as with the cmh-test")
(options, args) = parser.parse_args()

baselist,derivedlist=parse_comparestring(options.comp)

for chr,pos,mami,s in SyncReaderMajMin(options.sync):
        baseentries,derivedentries=parse_line(s,baselist,derivedlist)
        assstat=get_associastat(baseentries,derivedentries)
        
        topr=[chr,str(pos),str(assstat)]
        print "\t".join(topr)
        