#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from syncIO import SyncReaderMajMin, SyncWindowReader




def parse_line(syncs,baselist,derivedlist):
	# (C,A),(108,45),(90,47))
	dera=[]
	basea=[]
	for chr,pos,mami,parsed in syncs:
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
		dera.append(derc)
		basea.append(basec)
	return (basea,dera)



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


def computeDivergenceSum(count1,count2):
	assert(len(count1)==len(count2))
	
	disum=0.0
	for c1,c2 in zip(count1,count2):
		disum+=computeDivergence(c1,c2)
	return disum

def computeHeterozygositySum(count):
	hsum=0.0
	for c in count:
		hsum+=computeHeterozygosity(c)
	return hsum
	



def get_hscwindow(basec,derc):
	hs=computeHeterozygositySum(derc)
	d=computeDivergenceSum(basec,derc)
	if d<0.00000000000001:
		return 0.0
	return hs/d

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

for chr,start,end,syncs in SyncWindowReader(SyncReaderMajMin(options.sync),1000):
	if len(syncs)==0:
		continue
	
	basec,derc=parse_line(syncs,baselist,derivedlist)
	hsc=get_hscwindow(basec,derc)
	
	topr=[chr,str(start),str(hsc)]
	print "\t".join(topr)
