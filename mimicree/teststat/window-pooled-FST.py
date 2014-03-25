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
		
		dercov=float(derc[0]+derc[1])
		basecov=float(basec[0]+basec[1])
		derf=[float(derc[0])/dercov,float(derc[1])/dercov]
		basef=[float(basec[0])/basecov,float(basec[1])/basecov]
		dera.append(derf)
		basea.append(basef)
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

def computeHeterozygosity(f):
	h=1-f[0]**2.0-f[1]**2.0
	return h

def computeHeterozygositySum(freq):
	hsum=0.0
	for f in freq:
		hsum+=computeHeterozygosity(f)
	return hsum
	
def get_averagefreq(basef,derf):
	assert(len(basef)==len(derf))
	averagefreq=[]
	for bf,df in zip(basef,derf):
		avma=(bf[0]+df[0])/2.0
		avmi=(bf[1]+df[1])/2.0
		averagefreq.append([avma,avmi])
	return averagefreq

def get_pooledfst(basef,derf):
	avf=get_averagefreq(basef,derf)
	htotal=computeHeterozygositySum(avf)
	h1=computeHeterozygositySum(basef)
	h2=computeHeterozygositySum(derf)
	hav=(h1+h2)/2.0
	fst=(htotal-hav)/htotal
	return fst


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
	
	basef,derf=parse_line(syncs,baselist,derivedlist)
	pooledfst=get_pooledfst(basef,derf)
	
	topr=[chr,str(start),str(pooledfst)]
	print "\t".join(topr)
