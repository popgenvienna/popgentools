#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from syncIO import SyncReaderMajMin



        

def parse_line(entries,baselist,derivedlist):
	# (C,A),(108,45),(90,47))
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

"""
def nover2(n):
	nf=float(n)
	ret=(nf*(nf-1.0))/2.0
	return ret

def computeHeterozygosit(maj,min):
	cov=maj+min
	h=(nover2(cov)-nover2(maj)-nover2(min))/nover2(cov)
	return h
"""

def computeHeterozygosity(c):
	cov=c[0]+c[1]
	if cov==0:
		return 0.0
	fcov=float(cov)
	fmaj=float(c[0])/fcov
	fmin=float(c[1])/fcov
	h=1-fmaj**2.0-fmin**2.0
	return h

def computeDivergence(c1,c2):
	cov1=c1[0]+c1[1]
	cov2=c2[0]+c2[1]
	if cov1==0 or cov2==0:
		return 0.0
	sumcomp=float(cov1*cov2)
	summaj=float(c1[0]*c2[0])
	summin=float(c1[1]*c2[1])
	div=(sumcomp-summaj-summin)/sumcomp
	return div


def get_hsc(basec,derc):
	hs=computeHeterozygosity(derc)
	c=computeDivergence(basec,derc)
	if c<0.00000000000001:
		return 0.0
	return hs/c

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
	# [2L,15,(C,A),(108,45),(90,47))
	basec,derc=parse_line(s,baselist,derivedlist)
	hsc=get_hsc(basec,derc)
	
	topr=[chr,str(pos),str(hsc)]
	print "\t".join(topr)
        