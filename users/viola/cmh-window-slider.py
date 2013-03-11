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
                pval=a[-1]
                #2L_sim_w_su1    177693  A       234:0:0:0:0:0   174:0:0:0:0:0   236:0:0:1:0:0   176:0:0:0:0:0   133:0:0:0:0:0   194:0:0:0:0:0   284:0:0:0:0:0   166:0:0:0:0:0   108:0:4:0:0:0   40:0:2:0:0:0    NaN
                if(pval=="NaN"):
                        return (None,None,None)
                pval=float(pval)
                if not self.__logpvalue:
                        if(pval==0):
                                pval=310 # max
                        else:
                                pval= -1.0 * math.log10(pval)
                return(chr,pos,pval)
                
def printwindow(chr,startpos,endpos,pvals):
        midpos = (startpos+endpos)/2.0
        count = len(pvals)
        if count==0:
                return None
        sum = math.fsum(pvals) # geometric mean is the sum of the logs divided by n (usually the n's root but taking the log results in division by n)
        geomean=sum/float(count)
        print "{0}\t{1}\t{2}\t{3}".format(chr,int(midpos),count,geomean)
        

# /Users/robertkofler/dev/testdata/CMH_output_with_NaN.cmh.txt
# --cmh-file /Volumes/Temp1/simcoal2/sim_haplotype_diversity/sim-res/sync/cmh-g60-hdall/N1000-hd00.cmh --window-size 100 

parser = OptionParser()
parser.add_option("--cmh-file", dest="cmh", help="A file containing the output of the cmh-test")
parser.add_option("--window-size",dest="windowsize", help="The windowsize for the sliding window")
parser.add_option("--logpvalue", dest="logpvalue", default=False,action="store_true", help="Is the pvalue of the cmh test in -log10?")
(options, args) = parser.parse_args()

logpvalue=options.logpvalue
winsize=int(options.windowsize)
parser=CMHParser(logpvalue)

activechr=None
startpos=0
endpos=winsize
pvalsforwindow=[]

for line in open(options.cmh):
        chr,pos,pval=parser.parseLine(line)
        if chr is None:
                # skip invalid entries, eg p-value=NaN
                continue
        
        if activechr is None:
                activechr=chr
                
        if activechr !=chr:
                # a new chromosome has been encountered
                printwindow(activechr,startpos,endpos,pvalsforwindow)
                activechr=chr
                startpos=0
                endpos=winsize
                pvalsforwindow=[]
        
        if pos < startpos:
                # check if the the file is sorted; This pos< startpos should not happen if it is sorted
                raise ValueError("cmh file needs to be sorted by start position for" +chr + " " +activechr+" "+str(pos)+" " +str(startpos))
                
        if(pos>endpos):
                # the position of the SNP exceeds the window position -> print and reshift window
                printwindow(activechr,startpos,endpos,pvalsforwindow)
                while(pos>endpos):
                        startpos+=winsize
                        endpos+=winsize
                pvalsforwindow=[]
                
        # in any case the new pvalue needs to be appended somewhere        
        pvalsforwindow.append(pval)
        
# Also do not forget to print the last window
printwindow(activechr,startpos,endpos,pvalsforwindow)
        
                
        
