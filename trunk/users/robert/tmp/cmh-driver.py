#!/usr/bin/env python
import os
import sys


replicates=[1,3,5,10,15,20]
generations=[10,20,30,40,50,60]
cmhpath=sys.argv[1]


def generate_cmhcommand(cmhpath,input,output,rep,gen):
    samples=(10,20,30,40,50,60)
    samplehash={}
    c=1
    for s in samples:
        samplehash[s]=c
        c+=1

    baseoffset=len(samples)+1
    sampleoffset=samplehash[gen]

    tocomp=[]
    for i in range(0,rep):
        base=1+baseoffset*i
        der=base+sampleoffset
        tocomp.append("{0}-{1}".format(base,der))
    if rep==1:
        tocomp.append(tocomp[0])

    temppr=",".join(tocomp)
    topr="perl {0} --min-count 1 --min-coverage 1 --max-coverage 100000 --min-logpvalue 0.0 --population {1} --input {2} --output {3}".format(cmhpath,temppr,input, output)
    return topr


for r in replicates:
    for g in generations:
        dirname="cmh_r%s_g%s"%(r,g)
        if not os.path.exists(dirname):
            os.mkdir(dirname)

commandlist=[]
for r in replicates:
    for g in generations:
        dirname="cmh_r%s_g%s"%(r,g)
        for i in range(1,11):
            input="n%s.sync" % i
            output="%s/n%s.cmh"%(dirname,i)
            cmhcommand=generate_cmhcommand(cmhpath,input,output,r,g)
            commandlist.append(cmhcommand)
"""        
for c in commandlist:
    print c
"""



def submit_job_max_len(commandlist, max_processes):
    import subprocess
    import time
    sleep_time = 3.0
    processes = list()
    for command in commandlist:
        print 'running {n} processes. Submitting {proc}.'.format(n=len(processes),proc=str(command))
        processes.append(subprocess.Popen(command, shell=True, stdout=None))
        while len(processes) >= max_processes:
            time.sleep(sleep_time)
            processes = [proc for proc in processes if proc.poll() is None]
    while len(processes) > 0:
        time.sleep(sleep_time)
        processes = [proc for proc in processes if proc.poll() is None]

submit_job_max_len(commandlist, max_processes=8)
print "Done"

