#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from Diploids import Diploid
from Population import PopGenerator, Population
from FitnessFunctions import FitnessFunctionParser


parser = OptionParser()
parser.add_option("--Ne", dest="ne", help="the number of diploid individuals")
parser.add_option("--selection",action="append", dest="selection",help="A selection regime starting at some timepoint; of the form '--selection generation:s1:s2:h1:h2' or with epistasis '--selection generation:s1:s2:h1:h2:e12:eh'; e12=epistatic interaction term between SNP1 and SNP2; eh=epistatic heterozygosity only 1 or 2 (1..AaBb is already epistatic; 2..only AABB is epistatic")
parser.add_option("--r",dest="r",help="recombination rate between the two linked loci")
parser.add_option("--p1",dest="p1", help="starting frequency of first SNP")
parser.add_option("--p2",dest="p2", help="starting frequency of second SNP")
parser.add_option("--rsquared",dest="rsquared", help="LD at start")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--max-generations",dest="maxgen")

(options, args) = parser.parse_args()
repsim  = int(options.repsim)
twone   = int(options.ne)*2
r       = float(options.r)
p1      = float(options.p1)
p2      = float(options.p2)
rsquared    = float(options.rsquared)
maxgen = float(options.maxgen)
assert eh==1 or eh==2

selectiondictionary=FitnessFunctionParser.get_fitnessFunctionDictionary(options.selection)

for i in range(0,repsim):
        sA="S"
        sB="S"
        sAB="S"
        genA=int(maxgen)
        genB=int(maxgen)
        genAB=int(maxgen)
        pop=PopGenerator.ini_ld(twone,p1,p2,rsquared)

        counter=0
        ff=None
        while(not pop.is_fixed()):
                if (counter+1) in selectiondictionary:
                        ff=selectiondictionary[(counter+1)]
                pop=pop.getNextGeneration(twone,ff,r)       
                counter+=1
                if(sA =="S" and pop.is_fixedA()):
                        sA=pop.status(pop.countA())
                        genA=counter
                if(sB =="S" and pop.is_fixedB()):
                        sB=pop.status(pop.countB())
                        genB=counter
                if(sAB=="S" and pop.is_fixedAB()):
                        sAB=pop.status(pop.countAB())
                        genAB=counter
                if counter>=maxgen:
                        break    
        print "A\t{0}\t{1}\t{2}".format(sA,genA,i+1)
        print "B\t{0}\t{1}\t{2}".format(sB,genB,i+1)
        print "AB\t{0}\t{1}\t{2}".format(sAB,genAB,i+1)



        


