#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math

class Diploid:
        def __init__(self,hap1,hap2):
                """
                hap=(0,0)
                hap=(1,1)
                hap=(1,0)
                hap=(0,1)
                First is the selected one, second is the second selected one
                """
                self.__hap1=hap1
                self.__hap2=hap2
        
        def countab(self):
                ab=0
                if(self.__hap1[0]==0 and self.__hap1[1]==0):
                        ab+=1
                if(self.__hap2[0]==0 and self.__hap2[1]==0):
                        ab+=1
                return ab
        
        def countAB(self):
                AB=0
                if(self.__hap1[0]==1 and self.__hap1[1]==1):
                        AB+=1
                if(self.__hap2[0]==1 and self.__hap2[1]==1):
                        AB+=1
                return AB
        
        def countA(self):
                return (self.__hap1[0]+self.__hap2[0])
                

        def countB(self):
                return (self.__hap1[1]+self.__hap2[1])
        
        def __get_indici(self,r):
                """
                Whole recombination is happening here
                """
                firstindex=True #1
                if random.random()<0.5:
                        firstindex=False #0
                
                secondindex=firstindex
                if(random.random()<r):
                        secondindex=not firstindex
                return (firstindex,secondindex)
                
        def get_genotypeA(self):
                return self.countA()
                
        def get_genotypeB(self):
                return self.countB()
                
        def get_haplotypes(self):
                toret=[]
                toret.append(self.__hap1)
                toret.append(self.__hap2)
                return toret
        
        def get_gamete(self,r):
                firstindex,secondindex=self.__get_indici(r)
                # firstindex = False or True
                # secondindex = False or True
                
                g1=0
                if(firstindex):
                        # Firstindex = True = 1
                        g1=self.__hap1[0]
                else:
                        # Firstindex = False = 0
                        g1=self.__hap2[0]
                
                g2=0
                if(secondindex):
                        g2=self.__hap1[1]
                else:
                        g2=self.__hap2[1]
                return (g1,g2)
 