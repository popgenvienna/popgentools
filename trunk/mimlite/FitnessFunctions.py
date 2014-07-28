#!/usr/bin/env python
# developed by Robert Kofler
import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math


        

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
