#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math
from Population import PopGenerator, Population
from FitnessFunctions import FitnessFunctionParser
                   

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
        plt.title("Trajectories of recessive haplotype (ab)")
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
parser.add_option("--r",dest="r",help="recombination rate between the two linked loci")
parser.add_option("--selection",action="append", dest="selection",help="A selection regime starting at some timepoint; of the form '--selection generation:s1:s2:h1:h2' or with epistasis '--selection generation:s1:s2:h1:h2:e12:eh'; e12=epistatic interaction term between SNP1 and SNP2; eh=epistatic heterozygosity only 1 or 2 (1..AaBb is already epistatic; 2..only AABB is epistatic")
parser.add_option("--p1",dest="p1", help="p1 at the start")
parser.add_option("--p2",dest="p2", help="p2 at the start")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--max-generations",dest="maxgen")
parser.add_option("--output",dest="output")

(options, args) = parser.parse_args()
repsim = int(options.repsim)
twone  = int(options.ne)*2
r      = float(options.r)
p1     = float(options.p1)
p2     = float(options.p2)
maxgen = float(options.maxgen)
output = options.output

assert(p1+p2<=1.0)

selectiondictionary=FitnessFunctionParser.get_fitnessFunctionDictionary(options.selection)

trajectories_selected1=[]
trajectories_selected2=[]
trajectories_ab=[]
ld_decay=[]

for i in range(0,repsim):
        pop=PopGenerator.ini_competition(twone,p1,p2)
        freqsA=[]
        freqsB=[]
        ld=[]
        freqsab=[]
        freqsA.append(pop.get_frequencyA())
        freqsB.append(pop.get_frequencyB())
        freqsab.append(pop.get_frequencyab())
        ld.append(pop.get_rsquared())
        counter=0
        ff=None
        while(not pop.is_fixed()):
                if (counter+1) in selectiondictionary:
                        ff=selectiondictionary[(counter+1)]
                pop=pop.getNextGeneration(twone,ff,r)
                freqsA.append(pop.get_frequencyA())
                freqsB.append(pop.get_frequencyB())
                freqsab.append(pop.get_frequencyab())
                ld.append(pop.get_rsquared())
                if counter>=maxgen:
                        break        
                counter+=1
        trajectories_selected1.append(freqsA)
        trajectories_selected2.append(freqsB)
        trajectories_ab.append(freqsab)
        ld_decay.append(ld)

pad_trajectories(trajectories_selected1)
pad_trajectories(trajectories_selected2)
pad_trajectories(ld_decay)
pad_trajectories(trajectories_ab)

visualize_trajectories(trajectories_selected1, trajectories_selected2,trajectories_ab, ld_decay, output)


        


