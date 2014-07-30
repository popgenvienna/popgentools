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

def visualize_phenotypes(phenotypes,outputfile):
        import matplotlib.pyplot as plt
        import math
        simrep=len(phenotypes)

        plotcount=simrep
        rows=math.ceil(float(plotcount)/2.0)
        tmp=str(int(rows))+"2"
        plt.figure(figsize=(20,int(rows)*7))
        for k in range(1,plotcount+1):
                p=phenotypes[k-1]
                pc=len(p)
                rp=reduce_generations(p,int(pc/15))
                sbarg=tmp+str(k)
                plt.subplot(sbarg)
                plt.ylim(0,1.0)
                plt.boxplot(rp)
        plt.savefig(outputfile, bbox_inches=0)

                
def visualize_trajectories(trajectories,output):
        import matplotlib.pyplot as plt
        import math

        plotcount=len(trajectories)
        rows=math.ceil(float(plotcount)/2.0)
        tmp=str(int(rows))+"2"
        plt.figure(figsize=(20,int(rows)*7))
        #112

        for k in range(1,plotcount+1):
                s=trajectories[k-1]
                xachsis=range(1,s.generations()+1)
                sbarg=tmp+str(k)
                plt.subplot(sbarg)
                plt.ylim(0,1.0)
                plt.title(s.snpinfo())
                for i in range(0,s.repeats()):
                        plt.plot(xachsis,s.trajectory(i))

        if output is None:
                plt.show()
        else:
                plt.savefig(output, bbox_inches=0)


parser = OptionParser()
parser.add_option("--Ne", dest="ne", help="the number of diploid individuals")
parser.add_option("--selection",dest="selection",help="specifying the selection in the form fitmin:fitmax:mean:stdDev")
parser.add_option("--variant-contributions",dest="varcontri",help="A comma separated list specifiying the contribution of every variant to the phenotype")
parser.add_option("-p",dest="p", help="a comma separated list of starting allele frequencies")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--max-generations",dest="maxgen")
parser.add_option("--output-traj",dest="outputtraj",default=None,help="The output file")
parser.add_option("--output-phenotypes",dest="outputpheno",default=None,help="Boxplots of the phenotypes")
(options, args) = parser.parse_args()
repsim = int(options.repsim)
twone = int(options.ne) * 2
p = map(float,options.p.split(","))                       # start allele frequencies
selection = map(float,options.selection.split(":"))        # selection coefficient, mean, standard deviation
phenocontri = map(float, options.varcontri.split(","))               # contribution of the variants to the phenotype
assert(len(phenocontri) == len(p))
snpcount=len(phenocontri)
maxgen = int(options.maxgen)
startc = [int(twone*i) for i in p]

if(options.outputpheno is not None and repsim>9):
        raise ValueError("Can not print phenotypes for more than 9 simulations; either do not plot phenotypes or use less simulations")

fitnesCalc=FitnessCalculator(selection[0],selection[1],selection[2],selection[3])

# initialize the trajectory array; I fucking hate 3D arrays...
trajectories=[]
for i in range(0,snpcount):
        trajectories.append(TrajectoriesForSNP("{0}:{1}".format(i+1,phenocontri[i]),repsim))

phenotypes=[]
for i in range(0,repsim):
        phenotypes.append([])

for rep in range(0,repsim):
        # initialize
        pop=PopGenerator.ini_complete_linkage(twone,startc,phenocontri)
        for i in range(0,snpcount):
                freq=pop.get_frequencyAt(i)
                trajectories[i].appendFreq(rep,freq)
        phenotypes[rep].append(pop.get_relative_phenotypes())

        for g in range(0, maxgen):
                pop=pop.getNextGeneration(twone,fitnesCalc,phenocontri)
                phenotypes[rep].append(pop.get_relative_phenotypes())
                for i in range(0,snpcount):
                        freq=pop.get_frequencyAt(i)
                        trajectories[i].appendFreq(rep,freq)
visualize_trajectories(trajectories,options.outputtraj)
if(options.outputpheno is not None):
        visualize_phenotypes(phenotypes,options.outputpheno)
  


