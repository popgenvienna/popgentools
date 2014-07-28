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


def visualize_trajectories(trajectories_selected,trajectories_neutral, trajectories_AB, ld,output):
        import matplotlib.pyplot as plt
        import math

        plt.figure(figsize=(20,14))
        plt.subplot("221")
        plt.ylim(0,1.0)
        plt.title("Trajectories of selected SNPs (A)")
        for i in range(0,len(trajectories_selected)):
                at=trajectories_selected[i] #active trajectory
                xachsis=range(1,len(at)+1)
                plt.plot(xachsis,at)
                
        plt.subplot("222")
        plt.ylim(0,1.0)
        plt.title("Trajectories of neutral linked SNPs (B)")
        for i in range(0,len(trajectories_neutral)):
                at=trajectories_neutral[i] #active trajectory
                xachsis=range(1,len(at)+1)
                plt.plot(xachsis,at)
                
        plt.subplot("223")
        plt.ylim(0,1.0)
        plt.title("Trajectories of haplotype (AB)")
        for i in range(0,len(trajectories_AB)):
                at=trajectories_AB[i] #active trajectory
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
parser.add_option("--s",dest="s",help="selection coefficient of the beneficial allele")
parser.add_option("--r",dest="r",help="recombination rate between the two linked loci")
parser.add_option("--h",dest="het",help="heterozygosity")
parser.add_option("--p1",dest="p1", help="p1 at the start")
parser.add_option("--p2",dest="p2", help="p2 at the start")
parser.add_option("--rsquared",dest="rsquared",help="The extend of LD at the start")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--max-generations",dest="maxgen")
parser.add_option("--output",dest="output")

(options, args) = parser.parse_args()
repsim = int(options.repsim)
twone  = int(options.ne)*2
s      = float(options.s)
r      = float(options.r)
h      = float(options.het)
p1     = float(options.p1)
p2     = float(options.p2)
rsq    = float(options.rsquared)
maxgen = float(options.maxgen)
output = options.output

ff=FitnessFunctionNormal(s,h,0.0,0.5)
trajectories_selected=[]
trajectories_neutral=[]
trajectories_AB=[]
ld_decay=[]

for i in range(0,repsim):
        pop=PopGenerator.ini_ld(twone,p1,p2,rsq)
        freqsA=[]
        freqsB=[]
        freqsAB=[]
        ld=[]
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
        trajectories_selected.append(freqsA)
        trajectories_neutral.append(freqsB)
        trajectories_AB.append(freqsAB)
        ld_decay.append(ld)

pad_trajectories(trajectories_selected)
pad_trajectories(trajectories_neutral)
pad_trajectories(trajectories_AB)
pad_trajectories(ld_decay)

visualize_trajectories(trajectories_selected, trajectories_neutral, trajectories_AB, ld_decay, output)


        


