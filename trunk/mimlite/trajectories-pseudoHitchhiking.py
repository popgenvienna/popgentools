#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math
from Population import PopGenerator, Population
from FitnessFunctions import FitnessFunctionNormal
                   

def pad_trajectories(trajectories):
        maxlen=0
        for traj in trajectories:
                if(len(traj)>maxlen):
                        maxlen=len(traj)
        
        for traj in trajectories:
                lastelement=traj[-1]
                difference=maxlen-len(traj)
                for i in range(0,difference):
                        traj.append(lastelement)
        if(len(traj)!=maxlen):
                raise ValueError("failed padding")


def visualize_trajectories(trajectories_selected1,trajectories_selected2, trajectories_ab, ld,output):
        import matplotlib.pyplot as plt
        import math

        plt.figure(figsize=(20,14))
        plt.subplot("221")
        plt.ylim(0,1.0)
        plt.title("Trajectories of first selected SNPs (A)")
        for i in range(0,len(trajectories_selected1)):
                at=trajectories_selected1[i] #active trajectory
                xachsis=range(1,len(at)+1)
                plt.plot(xachsis,at)
                
        plt.subplot("222")
        plt.ylim(0,1.0)
        plt.title("Trajectories of second selected SNPs (B)")
        for i in range(0,len(trajectories_selected2)):
                at=trajectories_selected2[i] #active trajectory
                xachsis=range(1,len(at)+1)
                plt.plot(xachsis,at)

        plt.subplot("223")
        plt.ylim(0,1.0)
        plt.title("Trajectories of  haplotype (AB)")
        for i in range(0,len(trajectories_ab)):
                at=trajectories_ab[i] #active trajectory
                xachsis=range(1,len(at)+1)
                plt.plot(xachsis,at)
                
        plt.subplot("224")
        plt.ylim(0,1.0)
        plt.title("LD (r_squared)")
        for i in range(0,len(ld)):
                at=ld[i] #active trajectory
                xachsis=range(1,len(at)+1)
                plt.plot(xachsis,at)

        if output is None:
                plt.show()
        else:
                plt.savefig(output, bbox_inches=0)


parser = OptionParser()
parser.add_option("--Ne", dest="ne", help="the number of diploid individuals")
parser.add_option("--s1",dest="s1",help="selection coefficient of the beneficial allele of first selected SNP")
parser.add_option("--s2",dest="s2",help="selection coefficient of the beneficial allele of second selected SNP")
parser.add_option("--r",dest="r",help="recombination rate between the two linked loci")
parser.add_option("--h1",dest="het1",help="heterozygosity of first selected SNP")
parser.add_option("--h2",dest="het2",help="heterozygosity of second selected SNP")
parser.add_option("--p1",dest="p1", help="p1 at the start")
parser.add_option("--p2",dest="p2", help="p2 at the start")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--max-generations",dest="maxgen")
parser.add_option("--output",dest="output")

(options, args) = parser.parse_args()
repsim = int(options.repsim)
twone  = int(options.ne)*2
s1      = float(options.s1)
s2      = float(options.s2)
r      = float(options.r)
h1      = float(options.het1)
h2      = float(options.het2)
p1     = float(options.p1)
p2     = float(options.p2)
maxgen = float(options.maxgen)
output = options.output

assert(p1+p2<=1.0)

ff=FitnessFunctionNormal(s1,h1,s2,h2)

trajectories_selected1=[]
trajectories_selected2=[]
trajectories_AB=[]
ld_decay=[]

for i in range(0,repsim):
        pop=PopGenerator.ini_subfrequency(twone,p1,p2)
        freqsA=[]
        freqsB=[]
        ld=[]
        freqsAB=[]
        freqsA.append(pop.get_frequencyA())
        freqsB.append(pop.get_frequencyB())
        freqsAB.append(pop.get_frequencyAB())
        ld.append(pop.get_rsquared())
        counter=0
        while(not pop.is_fixed()):
                pop=pop.getNextGeneration(twone,ff,r)
                freqsA.append(pop.get_frequencyA())
                freqsB.append(pop.get_frequencyB())
                freqsAB.append(pop.get_frequencyAB())
                ld.append(pop.get_rsquared())
                if counter>=maxgen:
                        break        
                counter+=1
        trajectories_selected1.append(freqsA)
        trajectories_selected2.append(freqsB)
        trajectories_AB.append(freqsAB)
        ld_decay.append(ld)

pad_trajectories(trajectories_selected1)
pad_trajectories(trajectories_selected2)
pad_trajectories(ld_decay)
pad_trajectories(trajectories_AB)

visualize_trajectories(trajectories_selected1, trajectories_selected2,trajectories_AB, ld_decay, output)


        


