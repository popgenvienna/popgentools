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


def reduce_generations(phenotypes,foldreduction):
        reduced=[]
        for i, p in enumerate(phenotypes):
                if foldreduction>0 and i%foldreduction==0:
                        reduced.append(p)
        return reduced



                



parser = OptionParser()
parser.add_option("--Ne", dest="ne", help="the number of diploid individuals")
parser.add_option("--selection",action="append",dest="selection",help="specifying the selection in the form generation:fitmin:fitmax:mean:stdDev")
parser.add_option("--variant-contributions",dest="varcontri",help="A comma separated list specifiying the contribution of every variant to the phenotype")
parser.add_option("-p",dest="p", help="a comma separated list of starting allele frequencies")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--max-generations",dest="maxgen")

(options, args) = parser.parse_args()
repsim = int(options.repsim)
twone = int(options.ne) * 2
p = map(float,options.p.split(","))                       # start allele frequencies
phenocontri = map(float, options.varcontri.split(","))               # contribution of the variants to the phenotype
assert(len(phenocontri) == len(p))
snpcount=len(phenocontri)
maxgen = int(options.maxgen)
startc = [int(twone*i) for i in p]

selectiondictionary=FitnessCalculatorParser.get_fitnessFunctionDictionary(options.selection)




for rep in range(0,repsim):
        # initialize
        sN=["S" for i in range(0,snpcount)]
        genN=[maxgen for i in range(0,snpcount)]
        pop=PopGenerator.ini_complete_linkage(twone,startc,phenocontri)
        counter=0
        fitnesCalc=None
        while(not pop.is_fixed()):
                if (counter+1) in selectiondictionary:
                        fitnesCalc=selectiondictionary[(counter+1)]
                pop=pop.getNextGeneration(twone,fitnesCalc,phenocontri)
                counter+=1
                for i in range(0,snpcount):
                        if(sN[i] =="S" and pop.is_fixedAt(i)):
                                sN[i]=pop.get_statusAt(i)
                                genN[i]=counter
                if counter>=maxgen:
                        break

        for i in range(0,pop.countSNPs()):
                print "{0}\t{1}\t{2}\t{3}".format(i+1,sN[i],genN[i],rep+1)
                
