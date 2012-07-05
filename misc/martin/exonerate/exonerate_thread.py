import re
import os 
import collections
import sys 
from optparse import OptionParser, OptionGroup
import thread
import threading
import time
import multiprocessing

#Author: Martin Kapun
#version 2.0

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python exonerate_thread.py --path /output --fasta longest_isoform.fasta --query genome.fasta --par \"model p2g --refine region \" --maxthreads 12
2)	This script created two directories in the output folder: /fasta and /exonerate. The program splits the multifasta input in single fasta files (e.g. transcripts or proteins) and stores them in /fasta. The it starts multiple exonerate commands (one per single fasta). The number of parallel processes i defined with the parameter --maxthreads. To run the exonerate command, one has to define the query sequence and the parameters. The parameters need to be flanked by quote symbols (\"\") to avoid conflict with the parameters of the python script (e.g. \"--model p2g --bestn 1 --showalignment \"). The finished alignemnts will be saved to /exonerate. Optionally in the future, this script wwill produce a log file telling if the alignemnt worked.  
	""")

parser.add_option("-p", "--path", dest="path",help="Path to create /fast and /exonerate")
parser.add_option("-f", "--fasta", dest="fasta",help="multifasta")
parser.add_option("-m", "--maxthreads", dest="maxthreads", help="maximum number of parallel processes")
parser.add_option("-q", "--query", dest="query", help="query fasta")
parser.add_option("-a", "--par", dest="par", help="exonerate parameters")
parser.add_option("-e", "--merge", dest="merge", help="merge individual exonerate files and produce output? (\"y\" or \"n\") put \"o\" if you want output only")
parser.add_option("-s", "--split", dest="split", help="split multifasta in single files and create directory? (\"y\" or \"n\")")
parser.add_option("-l", "--list", dest="list", help="list of proteinnames or \"na\"")

parser.add_option_group(group)
(options, args) = parser.parse_args()	
	
#########################################################   CODE   #########################################################################
filenames=[]
if str(options.split)=="y" and str(options.merge)!="o":
	os.system("mkdir "+str(options.path)+"/fasta")
	os.system("mkdir "+str(options.path)+"/exonerate")
	for l in open(str(options.fasta),"r"):
		if ">" in l:
			outfile=str(options.path)+"/fasta/"+l.split()[0][1:]+".fasta"
			out=open(outfile,"w")
			out.write(l)
			filenames.append(l.split()[0][1:])
		if ">" not in l:
			out.write(l)
	def exonerate(file1, parameters):
		command="exonerate "+str(options.path)+"/fasta/"+file1+".fasta "+str(options.query)+" "+parameters+" >"+" "+str(options.path)+"/exonerate/"+file1+".2ex"
		os.system(command)
	def exonerate1(file1, parameters):
		command="exonerate "+str(options.path)+"/fasta/"+file1+".fasta "+str(options.query)+" "+parameters+" >"+" "+str(options.path)+"/exonerate/"+file1+".2ex"
		return command
	c=0	
	for item in filenames:
		threadcount=len(multiprocessing.active_children())
		while(threadcount > int(options.maxthreads)):
			print "Sleeping 5s Zzzz Zzzz Zzzzz"
			time.sleep(5)
			threadcount=len(multiprocessing.active_children())
		t=multiprocessing.Process(target=exonerate,args=(item, str(options.par).strip("\"")))
		c+=1
		print "executing: "+exonerate1(item, str(options.par).strip("\""))+" thread: "+str(c)
		t.start()
		
if str(options.split)=="n" and str(options.merge)!="o":
	os.system("mkdir "+str(options.path)+"/exonerate")
	os.system("mkdir "+str(options.path)+"/exonerate")
	if str(options.list)=="na":
		for subdir, dirs, files in os.walk(str(options.path)+"/fasta"):
			for a in files:
				filenames.append(a[:-6])
	else: 
		for l in open(options.list):
			filenames.append(l.rstrip())
	def exonerate(file1, parameters):
		command="exonerate "+str(options.path)+"/fasta/"+file1+".fasta "+str(options.query)+" "+parameters+" >"+" "+str(options.path)+"/exonerate/"+file1+".2ex"
		os.system(command)
	def exonerate1(file1, parameters):
		command="exonerate "+str(options.path)+"/fasta/"+file1+".fasta "+str(options.query)+" "+parameters+" >"+" "+str(options.path)+"/exonerate/"+file1+".2ex"
		return command
	c=0	
	for item in filenames:
		threadcount=len(multiprocessing.active_children())
		while(threadcount > int(options.maxthreads)):
			print "Sleeping 5s Zzzz Zzzz Zzzzz"
			time.sleep(5)
			threadcount=len(multiprocessing.active_children())
		t=multiprocessing.Process(target=exonerate,args=(item, str(options.par).strip("\"")))
		c+=1
		print "executing: "+exonerate1(item, str(options.par).strip("\""))+" thread: "+str(c)
		t.start()
	
	
if str(options.merge)=="y" or str(options.merge)=="o":
	if filenames==[]:
		for subdir, dirs, files in os.walk(str(options.path)+"/fasta"):
			for a in files:
				filenames.append(a[:-6])
	outfile2=str(options.path)+"/"+str(options.fasta).split("/")[-1]+".2exonerate"
	out2=open(outfile2,"w")
	outfile3=str(options.path)+"/2error.log"
	out3=open(outfile3,"w")
	outfile4=str(options.path)+"/proteins_error.2list"
	out4=open(outfile4,"w")
	test=0
	for item in filenames:
		try:
			open(str(options.path)+"/exonerate/"+item+".2ex","r")
		except:
			out3.write(str(options.path)+"/exonerate/"+item+".2ex is not present\n")	
			out4.write(item+"\n")
		else:
			for l in open(str(options.path)+"/exonerate/"+item+".2ex","r"):
				if "C4 Alignment:" in l:
					test=1
			if test==1:
				for l in open(str(options.path)+"/exonerate/"+item+".2ex","r"):
					out2.write(l)
				test=0
			else: 
				out3.write(str(options.path)+"/exonerate/"+item+".2ex is empty\n")
				out4.write(item+"\n")	
			
