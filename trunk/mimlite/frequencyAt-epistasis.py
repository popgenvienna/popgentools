#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from Diploids import Diploid
from Population import PopGenerator, Population
from FitnessFunctions import FitnessFunctionEpistasis


parser = OptionParser()
parser.add_option("--Ne", dest="ne", help="the number of diploid individuals")
parser.add_option("--s1",dest="s1",help="selection coefficient of the beneficial allele of first selected SNP")
parser.add_option("--s2",dest="s2",help="selection coefficient of the beneficial allele of second selected SNP")
parser.add_option("--r",dest="r",help="recombination rate between the two linked loci")
parser.add_option("--h1",dest="het1",help="heterozygosity of first selected SNP")
parser.add_option("--h2",dest="het2",help="heterozygosity of bad hitchhiking SNP")
parser.add_option("--eh",dest="eh",help="epistatic heterozygosity only 1 or 2 (1..AaBb is already epistatic; 2..only AABB is epistatic")
parser.add_option("--e12",dest="e12",help="epistatic interaction term between SNP1 and SNP2")
parser.add_option("--p1",dest="p1", help="starting frequency of first SNP")
parser.add_option("--p2",dest="p2", help="starting frequency of second SNP")
parser.add_option("--rsquared",dest="rsquared", help="LD at start")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--snapshots", dest="snapshots",default=None,help="A comma separated list of generations at which the frequency of the SNPs should be recorded")
(options, args) = parser.parse_args()
repsim  = int(options.repsim)
twone   = int(options.ne)*2
s1      = float(options.s1)
s2      = float(options.s2)
r       = float(options.r)
h1      = float(options.het1)
h2      = float(options.het2)
e12     = float(options.e12)
eh      = int(options.eh)
p1      = float(options.p1)
p2      = float(options.p2)
rsquared    = float(options.rsquared)
snapshots=set(map(int, options.snapshots.split(",")))
maxgen = max(snapshots)
assert eh==1 or eh==2

ff=FitnessFunctionEpistasis(s1,h1,s2,h2,e12,eh)
print  "# snp_id\tgeneration\tfrequency\treplicate"

for i in range(0,repsim):
        pop=PopGenerator.ini_ld(twone,p1,p2,rsquared)
        print "{0}\t{1}\t{2}\t{3}".format("A",0,pop.get_frequencyA() ,i+1)
        print "{0}\t{1}\t{2}\t{3}".format("B",0,pop.get_frequencyB() ,i+1)
        counter=0
        while(True):
                pop=pop.getNextGeneration(twone,ff,r)       
                counter+=1
                if counter in snapshots:
                        print "{0}\t{1}\t{2}\t{3}".format("A",counter,pop.get_frequencyA() ,i+1)
                        print "{0}\t{1}\t{2}\t{3}".format("B",counter,pop.get_frequencyB() ,i+1)
                if counter>=maxgen:
                        break   
 




        


