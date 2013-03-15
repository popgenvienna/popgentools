import sys
import random
import collections
import os
from optparse import OptionParser, OptionGroup
import re

def get_sensitivity(cTruePositives,cPositives):
	sens=float(cTruePositives)/float(cPositives)
	return sens


def get_specificity(cTrueNegatives,cNegatives):
	spec=float(cTrueNegatives)/float(cNegatives)
	return spec

def get_snpcount(cmhfile):
	counter=0
	for line in open(cmhfile):
		counter+=1
	return counter
	

def get_key(chr,pos):
	return str(chr)+":"+str(pos)

def get_selected_snps(file):
	posset=set({})
	for line in open(file):
		line=line.rstrip()
		a=line.split("\t")
		key=get_key(a[0],a[1])
		posset.add(key)
	return posset

parser = OptionParser()
parser.add_option("--cmh-file", 	dest="cmhfile",help="The cmh test file")
parser.add_option("--additive-file",dest="additivefile",help="The number of SNPs which should be randomly drawn")
(options, args) = parser.parse_args()

posset=get_selected_snps(options.additivefile)
cPositives=len(posset)
snpCount=get_snpcount(options.cmhfile)
cNegatives=snpCount-cPositives

cTested=0
cTruePositives=0


for line in open(options.cmhfile):
	line=line.rstrip()
	a=line.split("\t")
	key=get_key(a[0],a[1])
	cTested+=1
	if(key in posset):
		cTruePositives+=1
		
	fN=cPositives-cTruePositives
	tN=snpCount-cTested-fN
	
	sens=get_sensitivity(cTruePositives,cPositives)
	spec=get_specificity(tN,cNegatives)
	print "{0}\t{1}\t{2}\t{3}".format(cTested,cTruePositives,sens,spec)
	