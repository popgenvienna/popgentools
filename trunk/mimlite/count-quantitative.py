#!/usr/bin/env python
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
from scipy.stats import norm


class FitnessCalculator:
        def __init__(self,selCoef,mean,std):
                self.__selCoef=selCoef
                self.__mean=mean
                self.__std=std
                
                scale=norm.pdf(mean,loc=mean,scale=std)
                self.__scale=scale

                
        def getFitness(self,relativeFitness):
                sc=self.__selCoef
                scale=self.__scale
                pd=norm.pdf(relativeFitness,loc=self.__mean,scale=self.__std)
                
                fitness=pd*self.__selCoef/scale
                fitness+=1.0
                return fitness
                


class Diploid:
        def __init__(self,hap1,hap2,phenotypicContribution):
                """
                
                """
                assert(len(hap1)==len(hap2))
                assert(len(phenotypicContribution)==len(hap1))
                self.__hap1=hap1
                self.__hap2=hap2
                self.__phenoContri=phenotypicContribution
                self.__phenotype=self.__calcPhenotypicValue()
                self.__relativePhenotype=self.__calcRelativePhenotypicValue()

                
        def countSNPs(self):
                """
                get the SNP count
                """
                return len(self.__hap1)
        
        def genotypeAt(self,at):
                return (self.__hap1[at]+self.__hap2[at])

        def getPhenotypicValue(self):
                return self.__phenotype

        def getRelativePhenotypicValue(self):
                return self.__relativePhenotype
        
        def __calcPhenotypicValue(self):
                pheVal=0.0
                for i in range(0,self.countSNPs()):
                        genotAt=self.genotypeAt(i)
                        pc=self.__phenoContri[i]
                        if(genotAt==2):
                                pheVal+=pc
                        elif(genotAt==1):
                                pheVal+=float(pc)/2.0
                return pheVal

        def __calcRelativePhenotypicValue(self):
                maxPhenotype=0.0
                for pc in self.__phenoContri:
                        maxPhenotype+=pc
                iis=self.__phenotype
                relative=iis/maxPhenotype
                return relative
        
        def get_gamete(self):
                
                gamete=[]
                for i in range(0,self.countSNPs()):
                        r=random.random()
                        if(r<0.5):
                                gamete.append(self.__hap1[i])
                        else:
                                gamete.append(self.__hap2[i])
                return gamete
        

class Population:
        def __init__(self,diploids):
                self.__diploids=diploids
                self.__twone=float(2.0 * len(diploids))
                if(self.__twone==0):
                        raise ValueError("Invalid population size; must be larger than zero")
        def countSNPs(self):
                return self.__diploids[0].countSNPs()
                
        def get_countAt(self,n):
                sum=0
                for dipl in self.__diploids:
                        sum+=dipl.genotypeAt(n)
                return sum
       
        def get_frequencyAt(self,n):
                sum=self.get_countAt(n)
                return float(sum)/float(self.__twone)
                
        
        def get_statusAt(self,n):
                sum=self.get_countAt(n)
                twone=int(self.__twone)
                if(sum==0):
                        return "L"
                elif(sum==twone):
                        return "F"
                else:
                        return "S"
        
        def is_fixedAt(self,i):
                countn=self.get_countAt(i)
                if(countn>0 and countn<self.__twone):
                        return False
                else:
                        return True
                
        def is_fixed(self):
                twone=int(self.__twone)
                fixed=True
                for i in range(0,self.countSNPs()):
                        fixedAt=self.is_fixedAt(i)
                        if(not fixedAt):
                                fixed=False
                return fixed

        def get_relative_phenotypes(self):
                toret=[]
                for dip in self.__diploids:
                        toret.append(dip.getRelativePhenotypicValue())
                return toret
        
        def __binarysearch(self,fitnestuples,rand):
                hi=len(fitnestuples)
                lo=0
                mid=0
                while lo<hi:
                        mid=int((lo+hi)/2)
                        midval=fitnestuples[mid][0]
                        if(midval<rand):
                                lo=mid+1
                        elif(midval>rand):
                                hi=mid
                        else:
                                # if due to crazy chance the exact fitness value of a individual is hit, return the value of the next!
                                # zero belongs to the first thus the exact value belongs to the next
                                return fitnestuples[mid+1][1]
                return fitnestuples[lo][1]
        
        
        def getNextGeneration(self,twone,fitnesCalc,phenocontri):
                
                fitnestuples=[]
                fitsum=0
                for dipl in self.__diploids:
                        fit=fitnesCalc.getFitness(dipl.getRelativePhenotypicValue())
                        fitsum+=fit
                        fitnestuples.append((fitsum,dipl))
                
                ne=int(twone/2.0)
                novelpop=[]
                for i in range(0,ne):
                        r1,r2=(random.random()*fitsum,random.random()*fitsum)
                        dipl1,dipl2=self.__binarysearch(fitnestuples,r1),self.__binarysearch(fitnestuples,r2)
                        
                        gamete1=dipl1.get_gamete()
                        gamete2=dipl2.get_gamete()
                        novelpop.append(Diploid(gamete1,gamete2,phenocontri))
                return Population(novelpop)
    
class PopGenerator:
        
        @classmethod
        def ini_complete_linkage(cls,twone,startc,phenocontri):
                
                gametes=[]
                
                for sc in startc:
                        oc=twone-sc
                        gametecol=[1 for i in range(0,sc)] + [0 for i in range(0,oc)]
                        gametes.append(gametecol)
                
                
                diplcount=int(twone/2)
                diploids=[]
                for i in range(0,diplcount):
                        gamete1=[]
                        gamete2=[]
                        for col in gametes:
                                n=col.pop(int(random.random()*len(col)))
                                gamete1.append(n)
                        for col in gametes:
                                n=col.pop(int(random.random()*len(col)))
                                gamete2.append(n)
                        dipl=Diploid(gamete1,gamete2,phenocontri)
                        diploids.append(dipl)
                return Population(diploids)
                
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
parser.add_option("--selection",dest="selection",help="specifying the selection in the form selcoefficient:mean:stdDev")
parser.add_option("--variant-contributions",dest="varcontri",help="A comma separated list specifiying the contribution of every variant to the phenotype")
parser.add_option("-p",dest="p", help="a comma separated list of starting allele frequencies")
parser.add_option("--repeat-simulations",dest="repsim")
parser.add_option("--max-generations",dest="maxgen")

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

fitnesCalc=FitnessCalculator(selection[0],selection[1],selection[2])




for rep in range(0,repsim):
        # initialize
        sN=["S" for i in range(0,snpcount)]
        genN=[maxgen for i in range(0,snpcount)]
        pop=PopGenerator.ini_complete_linkage(twone,startc,phenocontri)
        counter=0
        while(not pop.is_fixed()):
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
                
