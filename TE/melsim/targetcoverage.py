import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math
import random


class TEinsert:
        def __init__(self,chr,pos,fwdrev,fam,popfreq,order,coverage,presencereads,overlap):
                self.chr=chr
                self.pos=pos
                self.fwdrev=fwdrev
                self.fam=fam
                self.popfreq=popfreq
                self.order=order
                self.coverage=coverage
                self.presencereads=presencereads
                self.overlap=overlap
                
def format_te(te):
        topr=[]
        topr.append(te.chr) #1
        topr.append(str(te.pos))
        topr.append(te.fwdrev)
        topr.append(te.fam)
        topr.append(str(te.popfreq))
        topr.append(te.order)
        topr.append("-")
        topr.append("-") # comment 8
        topr.append("100")
        topr.append("101")
        topr.append(str(te.popfreq))
        topr.append(str(te.coverage)) #12
        topr.append(str(te.presencereads)) #13
        topr.append(str(te.coverage-te.presencereads)) #14
        topr.append(str(te.overlap)) #15
        topr.append("-") #16 start range of rev insertion
        topr.append("-")
        topr.append("-")
        topr.append("-")
        topr.append("-")
        topr.append("-")
        topr.append("-")
        return "\t".join(topr)
        

def bootstrapreads(presence,cov,targetcov):
        presfrac=float(presence)/float(cov)
        
        count=0
        for i in range(0,targetcov):
                r=random.random()
                if r<presfrac:
                        count+=1
        return count
                
def load_tes(file):
        """
 1                                                                                      8                               12      13                                                               19     20
2L	905	FR	Quasimodo	0.994623655913978	LTR	-	ncorrdist=166	496	595	1	45	45	0	0	1215	1314	0.989247311827957	93	92	1	0
2L	1396	FR	Idefix	0.950774612693653	LTR	-	ncorrdist=224	1104	1203	0.988505747126437	87	86	1	0	1589	1688	0.91304347826087	46	42	4	0
2L	4866	FR	INE-1	0.974425287356322	TIR	-	ncorrdist=172	4564	4663	0.96551724137931	58	56	2	0	5069	5168	0.983333333333333	60	59	1	0
2L	13099.5	FR	INE-1	0.902352941176471	TIR	-	ncorrdist=199	12825	12924	0.84	50	42	8	0	13275	13374	0.964705882352941	85	82	3	0
Col9: start of the range of the forward insertion
Col10: end of the range of the forward insertion
Col11: population frequency estimated by the forward insertion
Col12: coverage of the forward insertion
Col13: TE-presence reads of the forward insertion
Col14: TE-absence reads of the reverse insertion
Col15: is the range of the forward insertion overlapping with a forward-range of another TE insertion (0..no; 1..yes)
        """
        tes=[]
        for l in open(file):
                l=l.rstrip("\n")
                a=l.split("\t")
                
                cov=0
                if a[11] !="-":
                        cov+=int(a[11])
                if a[18] !="-":
                        cov+=int(a[18])
                        
                presence=0
                if a[12] !="-":
                        presence+=int(a[12])
                if a[19] !="-":
                        presence+=int(a[19])
                overlap=0
                if a[14]=="1":
                        overlap=1
                if a[21]=="1":
                        overlap=1
                popfreq=float(presence)/float(cov)
                #  0    1       2       3       4     5           
                # (chr, pos, fwdrev, fam, popfreq,order,coverage,presencereads):
                te=TEinsert(a[0],a[1],a[2],a[3],popfreq,a[5],cov,presence,overlap)
                tes.append(te)
        return tes


parser = OptionParser()
parser.add_option("--tes",dest="tes",help="the basis genome; sample coverages from here")
parser.add_option("--target-cov",dest="targetcov",help="adjust coverages of this genome")
parser.add_option("--min-count",dest="mincount",help="the minimum count")
(options, args) = parser.parse_args()

mincount=int(options.mincount)
tes=load_tes(options.tes)
targetcov=int(options.targetcov)

for te in  tes:
        cov=te.coverage
        newcount=bootstrapreads(te.presencereads,cov,targetcov)
        if newcount < mincount:
                continue
        
        te.coverage=targetcov
        te.presencereads=newcount
        te.popfreq=float(newcount)/float(targetcov)
        assert(te.presencereads>0)
        print format_te(te)


        
