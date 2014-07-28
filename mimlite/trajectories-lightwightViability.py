#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections


class SelectionViability:
        def __init__(self,s,h):
                self.s=s
                self.h=h
                maxfit=max(1.0, 1+h*s, 1+s)
                self.survivalprob={0: (1.0/maxfit) , 1: (1+h*s)/maxfit , 2: (1+s)/maxfit}
        
        def nextGeneration(self,twone,population):
                diploids=population.getIndividuals()
                
                survivors=[]
                for dipl in diploids:
                        # let's kill the weak ones
                        genot=dipl.get_genotype()
                        survivalprob = self.survivalprob[genot]
                        if random.random()< survivalprob:
                                survivors.append(dipl)
                ne=int(twone/2.0)
                novelpop=[]
                survivorcount=len(survivors)
                for i in range(0,ne):
                        i1,i2=(int(random.random()*survivorcount) , int(random.random()*survivorcount))
                        while(i1==i2):
                                i2=int(random.random()*survivorcount)
 
                        d1,d2=survivors[i1],survivors[i2]
                        genotype=d1.get_gamete()+d2.get_gamete()
                        novelpop.append(Diploid(genotype))
                return Population(novelpop)
                

class Diploid:
        def __init__(self,genotype):
                self.__genotype=genotype
                
        def get_genotype(self):
                """
                genotype code
                0=w11 = AA
                1=w12 = Aa
                2=w22 = aa
                thus genotype=gamete1+gamete2
                """
                return self.__genotype
                # return 0=w11, 1=w12, or 2=w22
        
        def get_gamete(self):
                if(self.__genotype==0):
                        return 0
                elif(self.__genotype==2):
                        return 1
                elif(self.__genotype==1):
                        if(random.random()<0.5):
                                return 0
                        else:
                                return 1
                else:
                        raise ValueError("Invalid genotype")

class Population:
        def __init__(self,diploids):
                self.__diploids=diploids
                self.__twone=float(2.0 * len(diploids))
                if(self.__twone==0):
                        raise ValueError("Invalid population size; must be larger than zero")
                
        def getIndividuals(self):
                return self.__diploids
                

        @classmethod
        def initialize_population(cls,twone,ns):
                count_one=ns
                count_zero=twone-count_one
                haps=[0 for i in range(0,count_zero)]+[1 for i in range(0,count_one)]
                diploids=[]
                while(len(haps)>0):
                        draw_one=haps.pop(int(random.random()*len(haps)))
                        draw_two=haps.pop(int(random.random()*len(haps)))
                        diploids.append(Diploid(draw_one+draw_two)) # genotype is just the sum of the two gametes 0=0+0 1=0+1 2=1+1
                return Population(diploids)
        
        def get_frequency(self):
                twone=self.__twone
                genosum=0.0
                for d in self.__diploids:
                        genosum+=float(d.get_genotype())
                return genosum/twone
        
        def getNextGeneration(self,selectionregime,twone):
                nextGen=selectionregime.nextGeneration(twone,self)
                return nextGen

def visualize_trajectories(trajectories,output):
        import matplotlib.pyplot as plt
        import math

        plt.figure(figsize=(10,7))
        plt.ylim(0,1.0)
        plt.title("Trajectories of selected SNPs")
        for i in range(0,len(trajectories)):
                at=trajectories[i] #active trajectory
                xachsis=range(1,len(at)+1)
                plt.plot(xachsis,at)

        if output is None:
                plt.show()
        else:
                plt.savefig(output, bbox_inches=0)


parser = OptionParser()
parser.add_option("--Ne", dest="ne", help="the number of diploid individuals")
parser.add_option("-e",dest="h",help="the heterozygosity")
parser.add_option("-s",dest="s",help="the selection coefficient")
parser.add_option("-p",dest="p", help="start allele frquency")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--generations",dest="gen")
parser.add_option("--output",dest="output",default=None)
(options, args) = parser.parse_args()
repsim = int(options.repsim)
twone  = int(options.ne)*2
s      = float(options.s)
h      = float(options.h)
p      = float(options.p)
gen = int(options.gen)
startc = int(twone*p)
output=options.output

selfunc=SelectionViability(s,h)

trajectories=[]
for i in range(0,repsim):
        # Repeats of the experiment
        pop=Population.initialize_population(twone,startc)
        freqs=[]
        freqs.append(pop.get_frequency())

        for i in range(0,gen):
                # number of specified generations
                pop=pop.getNextGeneration(selfunc,twone)
                freqs.append(pop.get_frequency())
        trajectories.append(freqs)

visualize_trajectories(trajectories,output)

