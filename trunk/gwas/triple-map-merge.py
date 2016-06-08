#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import collections


 # realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

 # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../../bin")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)


class Measure:
	def __init__(self,chr,pos,measure):
		self.chr=chr
		self.pos=pos
		self.measure=measure


class MeasureReader:
	"""
	A light-weight pileup reader
	
	"""
	# chrset,args.col,args.trim
	def __init__(self,file,col,trim,debug):
		self.__filename=file
		self.__filehandle=open(file,"r")
		self.__col=col
		self.__trim=trim
		self.__buffer=[]
		self.__debug=debug
		self.__activechr=None

	def __iter__(self):
		return self
	
	def close(self):
		self.__filehandle.close()
		
	def buffer(self,tobuffer):
	    self.__buffer.append(tobuffer)
	    
	
	def spool2chr(self,chr):
		passed=set([])
		while(True):
			try:
				test=self.next()
			except StopIteration:
				return False
			if(test.chr==chr):
				self.buffer(test)
				self.__activechr=chr
				return True
			elif(self.__debug):
				if test.chr not in passed:
					print "Spooling "+test.chr
					passed.add(test.chr)
		return False
	
	def get_pos(self,pos):
		while(True):
			try:
				test=self.next()
			except StopIteration:
				return None
			if(test.pos==pos and test.chr==self.__activechr):
				return test
			elif(test.pos>pos or test.chr != self.__activechr):
				self.buffer(test)
				return None

		return None


	def next(self):
		# first empty the buffer
		if len(self.__buffer)>0:
		    return self.__buffer.pop()
		
		# if empty read the next line    
		line=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":
				raise StopIteration
			line=line.rstrip('\n')
			if line.startswith("#"):
				continue
			if line != "":
				break
		
		a=line.split("\t")
		chr,pos,measure=a[0],a[1],a[self.__col]
		if(self.__trim!=""):
			measure=measure.lstrip(self.__trim)
		return Measure(chr,int(pos),float(measure))

def get_chrset(chr):
	if "," in chr:
		return set(chr.split(","))
	else:
		return set([chr])

def get_chrorderandmaxsizes(file,chrset):
	chrorder=[]
	mapped=set([])
	chrsize=collections.defaultdict(lambda:0)
	for l in open(file):
		a=l.rstrip("\n").split("\t")
		chr=a[0]
		if(chr not in mapped and chr in chrset):
			chrorder.append(chr)
			mapped.add(chr)
		if chr in chrset:
			pos=int(a[1])
			if(pos>chrsize[chr]):
				chrsize[chr]=pos
	return chrorder,chrsize

def _get_small(measures):
	toret=100000000000
	for m in measures:
		if m < toret:
			toret=m
	return toret
	

def _get_large(measures):
	toret=-100000000000
	for m in measures:
		if m>toret:
			toret=m
	return toret

def _get_m1(measures):
	return measures[0]

def _get_m2(measures):
	return measures[1]

def _get_m3(measures):
	return measures[2]
		
def get_measure(take):
	if(take=="small"):
		return _get_small
	elif(take=="large"):
		return _get_large
	elif(take=="map1"):
		return _get_m1
	elif(take=="map2"):
		return _get_m2
	elif(take=="map3"):
		return _get_m3
	else:
		raise Exception("Unknown methods")

parser = argparse.ArgumentParser(description="""           
Description
-----------
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Output

Authors
-------
    Robert Kofler
""")

parser.add_argument("--map1", type=str, required=True, dest="map1", default=None, help="Results of mapper 1")
parser.add_argument("--map2", type=str, required=True, dest="map2", default=None, help="Results of mapper 2")
parser.add_argument("--map3", type=str, required=True, dest="map3", default=None, help="Results of mapper 3")
parser.add_argument("--take", type=str, required=False, dest="take", default="small", help="the minimum coverage [map1 | map2 | map3 | small | large]")
parser.add_argument("--chromosomes", type=str, required=False, dest="chr", default="X,2L,2R,3L,3R,4", help="Report results for the given comma separated list of chromosomes")
parser.add_argument("--column", type=int, required=False, dest="col", default="-1", help="Which column should be used for merging? Counting starts at zero and -1 means last column, -2 the one before the last etc")
parser.add_argument("--trim", type=str, required=False, dest="trim", default="", help="Remove the given string at the start of chosen column, e.g.: '1:2=' may be provided when using PoPoolation2")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="The output file")
parser.add_argument("--debug", required=False, dest="debug", default=False, action="store_true", help="The output file")
args = parser.parse_args()

debug=args.debug
chrset=get_chrset(args.chr)

print "Reading first input file to determine order and size of the major chromosomes: "+args.map1
chrorder,chrmaxsize=get_chrorderandmaxsizes(args.map1,chrset)
print "Done, will use sort order of chromosomes "+str(chrorder)


take=args.take

mr1=MeasureReader(args.map1,args.col,args.trim,debug)
mr2=MeasureReader(args.map2,args.col,args.trim,debug)
mr3=MeasureReader(args.map3,args.col,args.trim,debug)

mc=get_measure(args.take)

ofh=open(args.output,"w")

print "Start parsing all three files "
for chr in chrorder:
	print "Merging chromosome "+ chr
	f1=mr1.spool2chr(chr)
	if(not f1):
		raise Exception("Did not find "+chr +" in file "+args.map1)
	f2=mr2.spool2chr(chr)
	if(not f2):
		raise Exception("Did not find "+chr +" in file "+args.map2)
	f3=mr3.spool2chr(chr)
	if(not f3):
		raise Exception("Did not find "+chr +" in file "+args.map3)
	
	for p in range(1,chrmaxsize[chr]+1):
		m1=mr1.get_pos(p)
		m2=mr2.get_pos(p)
		m3=mr3.get_pos(p)
		if(m1 is not None and m2 is not None and m3 is not None):
			m=mc([m1.measure,m2.measure,m3.measure])
			topr=[chr,str(p),str(m)]
			ofh.write("\t".join(topr)+"\n")
ofh.close()
			
