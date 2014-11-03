#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from scipy.stats import norm
from Quantitative import *


                
class TrajectoriesForSNP:
        def __init__(self,snpinfo,trajcount):
                trajectories=[]
                for i in range(0,trajcount):
                        trajectories.append([])
                self.__traj=trajectories
                self.__snpinfo=snpinfo
                self.__trajcount=trajcount
                
        def appendFreq(self,repeatnumber,frequency):
                self.__traj[repeatnumber].append(frequency)

        def generations(self):
                return len(self.__traj[0])

        def snpinfo(self):
                return self.__snpinfo

        def repeats(self):
                return self.__trajcount

        def trajectory(self,i):
                return self.__traj[i]




                



parser = OptionParser()
parser.add_option("--Ne", dest="ne", help="the number of diploid individuals")
parser.add_option("--selection",action="append",dest="selection",help="specifying the selection in the form generation:fitmin:fitmax:mean:stdDev")
parser.add_option("--variant-contributions",dest="varcontri",help="A comma separated list specifiying the contribution of every variant to the phenotype")
parser.add_option("-p",dest="p", help="a comma separated list of starting allele frequencies")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--snapshots", dest="snapshots",default=None,help="A comma separated list of generations at which the frequency of the SNPs should be recorded")
(options, args) = parser.parse_args()
repsim = int(options.repsim)
twone = int(options.ne) * 2
p = map(float,options.p.split(","))                       # start allele frequencies
phenocontri = map(float, options.varcontri.split(","))               # contribution of the variants to the phenotype
assert(len(phenocontri) == len(p))
snpcount=len(phenocontri)

startc = [int(twone*i) for i in p]
snapshots=set(map(int, options.snapshots.split(","))) 
maxgen = max(snapshots)

selectiondictionary=FitnessCalculatorParser.get_fitnessFunctionDictionary(options.selection)

# initialize the trajectory array; I fucking hate 3D arrays...

print  "# snp_id\tgeneration\tfrequency\treplicate"
for rep in range(0,repsim):
        # initialize
        pop=PopGenerator.ini_complete_linkage(twone,startc,phenocontri)
        for i in range(0,snpcount):
                freq=pop.get_frequencyAt(i)
                print "{0}\t{1}\t{2}\t{3}".format(i+1,0,freq,rep+1)
        
        fitnesCalc=None
        for g in range(0, maxgen):
                if (g+1) in selectiondictionary:
                        fitnesCalc=selectiondictionary[(g+1)]
                pop=pop.getNextGeneration(twone,fitnesCalc,phenocontri)
                if (g+1) in snapshots: 
                        for i in range(0,snpcount):
                                freq=pop.get_frequencyAt(i)
                                print "{0}\t{1}\t{2}\t{3}".format(i+1,g+1,freq,rep+1)

  


