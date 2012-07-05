#!/usr/bin/python
import sys
import os
import collections
import operator
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
metafor = importr('metafor')

class homogenity:
    
    @classmethod
    def pop_count(cls,input):
        fh=open(input,"r")
        line = fh.readline()
        col=line.split("\t")
        chro = col.pop(0)
        pos = col.pop(0)
        rc = col.pop(0)
        cmh_pvalue=col.pop(-1)
        popcount = len(col)
        return popcount
            
    @classmethod
    def read_populations(cls,population,input):
    
        pop_pair = population.split(",")
        pop = []
        for p in pop_pair:
            if p.find("-")==1:
                p1,p2 = p.split("-")
                p1=p1.strip()
                p2=p2.strip()
                if p1:
                    pop.append(p1)
                if p2:
                    pop.append(p2)
            else:
                sys.stderr.write("\nERROR --> population pair should give pop1-pop2 [Example: 1-3,5-7].\n\n")
                sys.exit(1)
    
        pop_sorted = sorted(pop)
        popcount = homogenity.pop_count(input)
        
        if int(pop_sorted[-1])>int(popcount):
            sys.stderr.write("\nERROR --> population count out of range. Only "+str(popcount)+" populations are available in input file, but population "+str(pop_sorted[-1])+" used \n\n")
            sys.exit(1)
        elif int(len(pop_sorted))<2:
            sys.stderr.write("\nERROR --> population count should be at least 2.\n\n")
            sys.exit(1)
            
        elif int(len(pop_sorted))%2!=0:
            sys.stderr.write("\nERROR --> populations should be given in pair [example: 1-4,2-5,3-6]\n\n")
            sys.exit(1)
        return pop_pair

    @classmethod
    def get_pop_list(cls,col):
        pop_list=[]
        bases = ["A","T","C","G","N","del"]
        for p in col:
            cov = p.split(":")
            pop_list.append(dict(zip(bases,cov)))
        
        return pop_list

    @classmethod
    def get_major_alleles(cls,pop_list):
    
        act=tct=cct=gct=nct=delct=0
        nuc_dict = {}  
        for pop in pop_list:
            act = int(act)+int(pop["A"])
            tct = int(tct)+int(pop["T"])
            cct = int(cct)+int(pop["C"])
            gct = int(gct)+int(pop["G"])
            nct = int(nct)+int(pop["N"])
            delct = int(delct)+int(pop["del"])
        
        nuc_dict["A"] = int(act)
        nuc_dict["T"] = int(tct)
        nuc_dict["C"] = int(cct)
        nuc_dict["G"] = int(gct)
        nuc_dict["N"] = int(nct)
        nuc_dict["del"] = int(delct)
    
        sorted_nuc_dict = sorted(nuc_dict.iteritems(), key=operator.itemgetter(1),reverse=True)
        major=sorted_nuc_dict[0][0]
        minor=sorted_nuc_dict[1][0]
        return (major,minor)


    @classmethod
    def woolf_test(cls,woolf_frequency,pop_pair):
        r=robjects.r
        robjects.r('''
        woolf = function(x) {
          
          x <- x + 1 / 2
          k <- dim(x)[3]
          or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
          w <-  apply(x, 3, function(x) 1 / sum(1 / x))
          1 - pchisq(sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2), k - 1)
        }
            ''')
           
        woolf = robjects.globalenv['woolf']
        #print(woolf.r_repr())
        
        #v = robjects.FloatVector([36.0, 140.0, 23.0, 2.0, 21.0, 130.0, 24.0, 3.0,11.0, 130.0, 10.0, 4.0])
        v = robjects.FloatVector(woolf_frequency)
        a = robjects.r['array'](v, dim = [2,2,int(len(pop_pair))])
        res = woolf(a)
        woolf_pvalue="na"
        woolf_pvalue =  list(res)[0]
        return woolf_pvalue


    @classmethod
    def BD_test(cls,BD_frequency):
       
        ai_el=[]
        bi_el=[]
        ci_el=[]
        di_el=[]
        
        for i in range(0,len(BD_frequency)):
            ai_el.append(BD_frequency[i][0])
            bi_el.append(BD_frequency[i][1])
            ci_el.append(BD_frequency[i][2])
            di_el.append(BD_frequency[i][3])
    
        r=robjects.r
        ai_el_r = robjects.FloatVector(ai_el)
        bi_el_r = robjects.FloatVector(bi_el)
        ci_el_r = robjects.FloatVector(ci_el)
        di_el_r = robjects.FloatVector(di_el)
        whatadd = robjects.FloatVector([0.50,0.50])
        whatto = robjects.StrVector(["all","all"])
        
        #res<- rma.mh(ai=for_bd[1,], bi=for_bd[2,], ci=for_bd[3,], di=for_bd[4,], add=c(1/2,1/2), to=c("all", "all"))
        res = r['rma.mh'](ai=ai_el_r, bi=bi_el_r, ci=ci_el_r, di=di_el_r,add=whatadd,to=whatto)
        #breslow-day test results
        bd= [res.rx2('TA')[0],res.rx2('TAp')[0]]    #test stat, pval for Breslow-Day test
        if bd[0] is robjects.NA_Logical:
            bd = ['na','na']
            
        # bdpval = res$BDp
        # return BD, BDp # bd[0] , bd[1]
        bd_pvalue="na"
        bd_pvalue = bd[1]
        return bd_pvalue

    @classmethod
    def SyncReader(cls,input,output,pop_pair):
        
        ofh = open(output,"w")
        
        
        fh=open(input,"r")
        for line in fh:
            line=line.rstrip().lstrip()
            if line:
                col=line.split("\t")
                chro = col.pop(0)
                pos = col.pop(0)
                rc = col.pop(0)
                cmh_pvalue=col.pop(-1)
                pop_list=[]
                
                pop_list = homogenity.get_pop_list(col)
                (major,minor) = homogenity.get_major_alleles(pop_list)

                woolf_frequency=[]
                BD_frequency=[]
                for p in pop_pair:
                    p1,p2 = p.split("-")
                    p1=int(p1.strip())-1
                    p2=int(p2.strip())-1
                    woolf_frequency.append(float(pop_list[p1][major]))
                    woolf_frequency.append(float(pop_list[p2][major]))
                    woolf_frequency.append(float(pop_list[p1][minor]))
                    woolf_frequency.append(float(pop_list[p2][minor]))
                    
                    BD_frequency.append([float(pop_list[p1][major]),float(pop_list[p1][minor]),float(pop_list[p2][major]),float(pop_list[p2][minor])])
    
                woolf_pvalue = homogenity.woolf_test(woolf_frequency,pop_pair)
                bd_pvalue = homogenity.BD_test(BD_frequency)
                
                pop_str = "\t".join(col)
                print >> ofh, "\t".join(map(str,[chro,pos,rc,pop_str,cmh_pvalue,woolf_pvalue,bd_pvalue]))
                #print "\t".join(map(str,[chro,pos,rc,pop_str,cmh_pvalue,woolf_pvalue,bd_pvalue]))
        ofh.close()
        fh.close()

