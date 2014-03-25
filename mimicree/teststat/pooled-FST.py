#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from syncIO import SyncReaderMajMin


        

def parse_line(entries,baselist,derivedlist):
	parsed=entries
	derc=[0,0]
	basec=[0,0]
	for i in derivedlist:
		key=i-1
		derc[0]+=parsed[key][0]
		derc[1]+=parsed[key][1]
	for i in baselist:
		key=i-1
		basec[0]+=parsed[key][0]
		basec[1]+=parsed[key][1]
	return (basec,derc)

def parse_comparestring(comp):
	a=comp.split(",")
	
	baselist=[]
	derivedlist=[]
	for e in a:
		b,d=map(int,e.split("-"))
		baselist.append(b)
		derivedlist.append(d)
	return (baselist,derivedlist)


def get_averagefreq(bf,df):
	avma=(bf[0]+df[0])/2.0
	avmi=(bf[1]+df[1])/2.0
	return (avma,avmi)

def computeHeterozygosity(f):
	h=1-f[0]**2.0-f[1]**2.0
	return h

def getFreq(c):
	cov=c[0]+c[1]
	if cov==0:
		return 0.0
	fcov=float(cov)
	fmaj=float(c[0])/fcov
	fmin=float(c[1])/fcov
	return (fmaj,fmin)

def get_pooled_fst(basec,derc):
	basef=getFreq(basec)
	derf=getFreq(derc)
	avf=get_averagefreq(basef,derf)
	htotal=computeHeterozygosity(avf)
	if htotal<0.0000000001:
		return 0.0
	h1=computeHeterozygosity(basef)
	h2=computeHeterozygosity(derf)
	hav=(h1+h2)/2.0
	fst=(htotal-hav)/htotal
	return fst



parser = OptionParser()
parser.add_option("--sync",dest="sync",help="A file containing the cmh results")
parser.add_option("--compare",dest="comp",help="A comparision string as with the cmh-test")
(options, args) = parser.parse_args()

baselist,derivedlist=parse_comparestring(options.comp)

for chr,pos,mami,s in SyncReaderMajMin(options.sync):
        basec,derc=parse_line(s,baselist,derivedlist)
        poolfst=get_pooled_fst(basec,derc)
        
        topr=[chr,str(pos),str(poolfst)]
        print "\t".join(topr)






        
        
                  
                  

  

#2L      6580    C       2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   2:0:70:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   2:0:70:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    0:0:118:0:0:0   0:0:72:0:0:0    2:0:116:0:0:0   0:0:72:0:0:0