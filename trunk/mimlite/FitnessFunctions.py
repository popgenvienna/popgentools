#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math

class FitnessFunctionParser:
        
        @classmethod
        def get_fitnessFunctionDictionary(cls,arguments):
                ffd={}
                for t in arguments:
                        a=t.split(":")
                           
                        generation=int(a[0])
                        s1=float(a[1])
                        s2=float(a[2])
                        h1=float(a[3])
                        h2=float(a[4])
                        if(len(a)==5):
                                ffd[generation]=FitnessFunctionNormal(s1,h1,s2,h2)
                        elif(len(a)==7):
                                e12=float(a[5])
                                eh=int(a[6])
                                assert eh==1 or eh==2
                                ffd[generation]=FitnessFunctionEpistasis(s1,h1,s2,h2,e12,eh)
                        else:
                                raise Error("Invalid arguments; fitness must either contain 5 fields (no epistasis) or 7 fields; "+a)
                if 1 not in ffd:
                        raise Error("Fitness starting at generation 1 needs to be provided; Example of neutral drift in the beginning '--selection 1;0;0;0;0' ")
                return ffd
                
        

class FitnessFunctionEpistasis:
        def __init__(self,s1,h1,s2,h2,e12,eh):
                self.__g2fA={0:1 , 1:1+h1*s1 , 2:1+s1}
                self.__g2fB={0:1 , 1:1+h2*s2 , 2:1+s2}
                self.__e12=e12
                self.__eh=eh
        def fitness(self,gA,gB):
                fitA=self.__g2fA[gA]
                fitB=self.__g2fB[gB]
                fit=fitA*fitB
                if(gA>=self.__eh and gB>=self.__eh):
                                fit=fit*(1+self.__e12)
                return fit
        
        
class FitnessFunctionNormal:
        def __init__(self,s1,h1,s2,h2):
                self.__g2fA={0:1 , 1:1+h1*s1 , 2:1+s1}
                self.__g2fB={0:1 , 1:1+h2*s2 , 2:1+s2}

        def fitness(self,gA,gB):
                fitA=self.__g2fA[gA]
                fitB=self.__g2fB[gB]
                fit=fitA*fitB
                return fit
