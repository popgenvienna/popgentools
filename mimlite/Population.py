#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math
from Diploids import Diploid

                
class Population:
        def __init__(self,diploids):
                self.__diploids=diploids
                self.__twone=float(2.0 * len(diploids))
                if(self.__twone==0):
                        raise ValueError("Invalid population size; must be larger than zero")
                countA=0
                countB=0
                for d in diploids:
                        countA+=d.countA()
                        countB+=d.countB()
                        
                self.__countA=countA
                self.__countB=countB
       
        def get_frequencyA(self):
                return float(self.__countA)/float(self.__twone)
        
        def get_frequencyB(self):
                return float(self.__countB)/float(self.__twone)
                
        def get_frequencyab(self):
                cab=self.countab()
                return float(cab)/float(self.__twone)
                
        def get_frequencyAB(self):
                cAB=self.countAB()
                return float(cAB)/float(self.__twone)
        
        def countAB(self):
                cAB=0
                for d in self.__diploids:
                        cAB+=d.countAB()
                return cAB
        
        def countA(self):
                return self.__countA
        
        def countB(self):
                return self.__countB

        def countab(self):
                cab=0
                for d in self.__diploids:
                        cab+=d.countab()
                return cab

        def status(self,count):
                twone=int(self.__twone)

                if(count==0):
                        return "L"
                elif(count==twone):
                        return "F"
                else:
                        return "S"
   
        def get_rsquared(self):
                """
                measure of ld
                """
                haplotypes=[]
                for dip in self.__diploids:
                        haplotypes.extend(dip.get_haplotypes())
                        
                p1,q1,p2,q2=(0,0,0,0)
                x11,x22,x12,x21=(0,0,0,0)
                
                for h in haplotypes:
                        if h ==(0,0):
                                x11+=1
                        elif h==(1,1):
                                x22+=1
                        elif h==(1,0):
                                x21+=1
                        elif h==(0,1):
                                x12+=1
                        else:
                                raise ValueError("Invalid haplotype")
                                
                        if(h[0]==1):
                                p1+=1
                        elif(h[0]==0):
                                q1+=1
                        else:
                                raise ValueError("Invalid allele frequency")
                        
                        if(h[1]==1):
                                p2+=1
                        elif(h[1]==0):
                                q2+=1
                        else:
                                raise ValueError("Invalid allele frequency")
                hapcount=len(haplotypes)
                p1,q1,p2,q2 = float(p1)/float(hapcount),float(q1)/float(hapcount),float(p2)/float(hapcount),float(q2)/float(hapcount)
                x11,x22,x12,x21 = float(x11)/float(hapcount),float(x22)/float(hapcount),float(x12)/float(hapcount),float(x21)/float(hapcount)
                D=x11*x22-x12*x21
                Dsquared=D**2
                div=p1*p2*q1*q2
                if(div==0):
                        return 0.0
                rsquared=Dsquared/div
                
                # stupid rounding mistakes...
                assert(rsquared<1.1 and rsquared>=0.0)
                
                if(rsquared>1.0):
                        rsquared=1.0
                if(rsquared<0.0):
                        rsquared=0.0
                return rsquared

        def is_fixedA(self):
                twone=int(self.__twone)
                counta=self.__countA
                if(counta>0 and counta<twone):
                        return False
                else:
                        return True

        def is_fixedB(self):
                twone=int(self.__twone)
                countb=self.__countB
                if(countb>0 and countb<twone):
                        return False
                else:
                        return True
        
        def is_fixedAB(self):
                twone=int(self.__twone)
                countAB=self.countAB()
                if(countAB>0 and countAB<twone):
                        return False
                else:
                        return True

        def is_fixedab(self):
                twone=int(self.__twone)
                countab=self.countab()
                if(countab>0 and countab<twone):
                        return False
                else:
                        return True
                
                
        def is_fixed(self):
                isfixedA=self.is_fixedA()
                isfixedB=self.is_fixedB()
                if(isfixedA and isfixedB):
                        return True
                else:
                        return False
        
        
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
        

        
        def getNextGeneration(self,twone,ff,r):
                fitnestuples=[]
                fitsum=0
                for dipl in self.__diploids:
                        fit=ff.fitness(dipl.get_genotypeA(),dipl.get_genotypeB())
                        fitsum+=fit
                        fitnestuples.append((fitsum,dipl))
                
                ne=int(twone/2.0)
                novelpop=[]
                for i in range(0,ne):
                        r1,r2=(random.random()*fitsum,random.random()*fitsum)
                        dipl1,dipl2=self.__binarysearch(fitnestuples,r1),self.__binarysearch(fitnestuples,r2)
                        
                        gamete1=dipl1.get_gamete(r)
                        gamete2=dipl2.get_gamete(r)
                        novelpop.append(Diploid(gamete1,gamete2))
                return Population(novelpop)
        
                        
class PopGenerator:

        @classmethod
        def ini_subfrequency(cls,twone,p1,p2):
                assert(p2<p1)
                
                sc=int(p2*twone)
                tc=int(p1*twone)-sc
                oc=twone-sc-tc
                gametecol=[(1,1) for i in range(0,sc)] +[(1,0) for i in range(0,tc)] + [(0,0) for i in range(0,oc)]
                diploids=[]
                while(len(gametecol)>0):
                        gamete1=gametecol.pop(int(random.random()*len(gametecol)))
                        gamete2=gametecol.pop(int(random.random()*len(gametecol)))
                        diploids.append(Diploid(gamete1,gamete2)) 
                return Population(diploids)
        
        @classmethod
        def ini_complete_linkage(cls,twone,p):

                sc=int(p*twone)
                oc=twone-sc
                gametecol=[(1,1) for i in range(0,sc)] + [(0,0) for i in range(0,oc)]
                
                diploids=[]
                while(len(gametecol)>0):
                        gamete1=gametecol.pop(int(random.random()*len(gametecol)))
                        gamete2=gametecol.pop(int(random.random()*len(gametecol)))
                        diploids.append(Diploid(gamete1,gamete2)) 
                return Population(diploids)

        @classmethod
        def ini_competition(cls,twone,p1,p2):
                assert(p1+p2<=1.0)
                if(p1+p2>1.0):
                        raise ValueError("Sum of starting allele frequencies need to be smaller than one")
                x12=int(twone*p1)
                x21=int(twone*p2)
                x11=twone-x12-x21
                        
                gametecol=[(0,0) for i in range(0,x11)] + [(0,1) for i in range(0,x21)] + [(1,0) for i in range(0,x12)]
                diploids=[]
                # x11 = (1,1) x22=(0,0) 
                while(len(gametecol)>0):
                        gamete1=gametecol.pop(int(random.random()*len(gametecol)))
                        gamete2=gametecol.pop(int(random.random()*len(gametecol)))
                        diploids.append(Diploid(gamete1,gamete2)) 
                return Population(diploids)

        
        @classmethod
        def ini_ld(cls,twone,p1,p2,rsquared):
               
                q1=1.0-p1
                q2=1.0-p2
                div=p1*p2*q1*q2
                D=math.sqrt(rsquared*div)
                x11=p1*p2+D
                x22=q1*q2+D
                x12=p1*q2-D
                x21=q1*p2-D
                
                x11,x22,x12=int(x11*twone),int(x22*twone),int(x12*twone)
                x21=twone-x11-x22-x12
                if(x11<0 or x22 < 0 or x12<0 or x21<0):
                        raise ValueError("Invalid value for rsquared; Impossible combination of allele frequencies and rsquared")
                

                gametecol=[(1,1) for i in range(0,x11)] + [(0,0) for i in range(0,x22)] + [(0,1) for i in range(0,x21)] + [(1,0) for i in range(0,x12)]
                diploids=[]
                # x11 = (1,1) x22=(0,0) 
                while(len(gametecol)>0):
                        gamete1=gametecol.pop(int(random.random()*len(gametecol)))
                        gamete2=gametecol.pop(int(random.random()*len(gametecol)))
                        diploids.append(Diploid(gamete1,gamete2)) 
                return Population(diploids)
                        




        


