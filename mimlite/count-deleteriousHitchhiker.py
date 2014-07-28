#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math
from Population import PopGenerator, Population
from FitnessFunctions import FitnessFunctionNormal

parser = OptionParser()
parser.add_option("--Ne", dest="ne", help="the number of diploid individuals")
parser.add_option("--s1",dest="s1",help="selection coefficient of the beneficial allele of first selected SNP")
parser.add_option("--s2",dest="s2",help="selection coefficient of the bad hitchhiking SNP (please provide a postive value)")
parser.add_option("--r",dest="r",help="recombination rate between the two linked loci")
parser.add_option("--h1",dest="het1",help="heterozygosity of first selected SNP")
parser.add_option("--h2",dest="het2",help="heterozygosity of bad hitchhiking SNP")
parser.add_option("--p",dest="p", help="starting frequency of SNPs")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--max-generations",dest="maxgen")

(options, args) = parser.parse_args()
repsim = int(options.repsim)
twone  = int(options.ne)*2
s1      = float(options.s1)
s2      = float(options.s2)
r      = float(options.r)
h1      = float(options.het1)
h2      = float(options.het2)
p     = float(options.p)
maxgen = float(options.maxgen)

ff=FitnessFunctionNormal(s1,h1,s2,h2)

for i in range(0,repsim):
        sA="S"
        sB="S"
        sAB="S"
        genA=int(maxgen)
        genB=int(maxgen)
        genAB=int(maxgen)
        pop=PopGenerator.ini_complete_linkage(twone,p)
        counter=0
        while(not pop.is_fixed()):
                pop=pop.getNextGeneration(twone,ff,r)    
                counter+=1
                
                if(sA =="S" and pop.is_fixedA()):
                        sA=pop.status(pop.countA())
                        genA=counter
                if(sB =="S" and pop.is_fixedB()):
                        sB=pop.status(pop.countB())
                        genB=counter
                if(sAB =="S" and pop.is_fixedAB()):
                        sAB=pop.status(pop.countAB())
                        genAB=counter
                if counter>=maxgen:
                        break    

        print "A\t{0}\t{1}\t{2}".format(sA,genA,i+1)
        print "B\t{0}\t{1}\t{2}".format(sB,genB,i+1)
        print "AB\t{0}\t{1}\t{2}".format(sAB,genAB,i+1)        




        


