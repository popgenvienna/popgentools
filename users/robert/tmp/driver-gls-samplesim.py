#!/usr/bin/env python
import os
import sys

# --min-count 5 --min-coverage 5 --max-coverage 100 --input sim-futsch2/Dmel-c30-91-92-93-105.sim.sync --output gls/Dmel-c30-91-92-93-105.mc5.simf2.gls.txt
path_gls="/Volumes/Temp2/Robert/popgentools/users/robert/pool-seq-paper/gls-pvalues.pl"
path_binomial="/Volumes/Temp2/Robert/popgentools/users/robert/sample-permute.py"
simulations=100

input=sys.argv[1]
outputdir=sys.argv[2]
outputdir=outputdir.rstrip("/")

commandpairs=[]
for i in range(1,simulations+1):
    outputfile="%s/temp%i"%(outputdir,i)
    samplecommand="python %s --input %s > %s"%s(path_binomial, input,outputfile)
    glsfile="%s/n%i.gls"%(outputdir,i)
    glscommand="perl %s --remove-temp --min-count 5 --min-coverage 5 --max-coverage 100 --input %s --output %s"%(path_gls,outputfile,glsfile)
    commandpairs.append([samplecommand,glscommand,i])
    

def run(args):
    run1=args[0]
    run2=args[1]
    counter=args[2]
    
    print "%i-1 Executung %s"%(counter,run1)
    os.system(run1)
    print "%i-2 Executing %s"%(counter,run2)
    os.system(run2)

pool = multiprocessing.Pool(processes = 8)
pool.map(run,commandpairs)
pool.close() # no more jobs accepted by pool
pool.join() # wait for jobs to finish

