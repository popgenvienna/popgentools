import sys
import collections
import math
import re
import random
from optparse import OptionParser, OptionGroup

class CMHParser:
        def __init__(self,logpvalue):
                self.__logpvalue=logpvalue
        
        def parseLine(self,line):
                l = line.rstrip()
                a=l.split("\t")
                chr=a[0]
                pos=int(a[1])
                pval=float(a[-1])
                if(pval=="NaN"):
                        return (None,None,None)
                
                if not self.__logpvalue:
                        if(pval==0):
                                pval=310 # max
                        else:
                                pval= -1.0* math.log10(pval)
                return(chr,pos,pval)
                
def printwindow(chr,startpos,endpos,pvals):
        midpos = (startpos+endpos)/2.0
        count = len(pvals)
        if count==0:
                return None
        sum = math.fsum(pvals) # geometric mean is the sum of the logs divided by n (usually the n's root but taking the log results in division by n)
        geomean=sum/float(count)
        print "{0}\t{1}\t{2}\t{3}".format(chr,int(midpos),count,geomean)
        
  

parser = OptionParser()
parser.add_option("--cmh-file", dest="cmh", help="A file containing the output of the cmh-test")
parser.add_option("--window-size",dest="windowsize", help="The windowsize for the sliding window")
parser.add_option("--logpvalue", dest="logpvalue", action="store_true", help="Is the pvalue of the cmh test in -log10?")
(options, args) = parser.parse_args()

winsize=int(options.windowsize)
parser=CMHParser(options.logpvalue)

activechr=None
startpos=0
endpos=winsize
pvalsforwindow=[]

for line in open(options.cmh):
        chr,pos,pval=parser.parseLine(line)
        if chr is None:
                continue
        
        if activechr is None:
                activechr=chr
                
        if activechr !=chr:
                printwindow(activechr,startpos,endpos,pvalsforwindow)
                activechr=chr
                startpos=0
                endpos=winsize
                pvalsforwindow=[]
        
        if pos < startpos:
                raise ValueError("cmh file needs to be sorted by start position for" +chr + " " +activechr+" "+str(pos)+" " +str(startpos))
                
        if(pos>endpos):
                printwindow(activechr,startpos,endpos,pvalsforwindow)
                while(pos>endpos):
                        startpos+=winsize
                        endpos+=winsize
                pvalsforwindow=[]
        # in any case the new pvalue needs to be appended somewhere        
        pvalsforwindow.append(pval)
printwindow(activechr,startpos,endpos,pvalsforwindow)
        
                
        
