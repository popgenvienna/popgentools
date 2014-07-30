#!/usr/bin/env python
from scipy.stats import norm
import random
import collections

class FitnessCalculator:
        def __init__(self,fitmin,fitmax,mean,std):
                self.__fitmin=fitmin
                self.__fitmax=fitmax
                self.__mean=mean
                self.__std=std
                scale=norm.pdf(mean,loc=mean,scale=std)
                self.__scale=scale

                
        def getFitness(self,relativeFitness):
                sc=self.__fitmax-self.__fitmin
                scale=self.__scale
                pd=norm.pdf(relativeFitness,loc=self.__mean,scale=self.__std)
                
                fitness=pd*sc/scale
                fitness+=self.__fitmin
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