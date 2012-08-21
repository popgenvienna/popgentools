import modules.RCMH
from modules.CMH import SyncReader,CMHReader
from optparse import OptionParser, OptionGroup
import collections
import rpy2.robjects as robjects
import math

#########################################################   HELP   #########################################################################
usage="python %prog --input input.sync --div output_cov100.div --pops 1,2,4 --min-count 2,2,3 --out output.af"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

based on a divergence dataset (--div) this script calculates the allele frequencies of the divergent allele from the species causing the contamination for each populations defined (--pops) in the input (--input). Additionally, you need to define a minimum count.  The output (--output) consists of columns for the chromosome, the position, the allele of the contaminating species followed by columns with allele frequencies for each population in (--pops). Per default the script is expecting a sync file as an input. However, it also handles CMH outputs (with a P-value at the end). Then, you need to put --type CMH

""") 

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--div", dest="div", help="divergence file")
parser.add_option("--pops", dest="pops", help="populations to test")
parser.add_option("--out", dest="out", help="output files")
parser.add_option("--mincount", dest="mincount", help="minimum count (for each population separately), e.g. 2,2,2")
parser.add_option("--type",dest="type",help="file type: cmh output or sync file",default="sync")

parser.add_option_group(group)
(options, args) = parser.parse_args()

########################################################### code ###################################################


################################################## read parameters #################################################


pops=map(int,options.pops.split(","))
out=open(options.out+"_table.txt","w")

############################################## test for format of input file ###################################

if options.type=="sync":
	filehandle=SyncReader(options.input)
else:
	filehandle=CMHReader(options.input)

############################### read divergence file #########################################################################

divhash=collections.defaultdict(lambda:"")
for l in open(options.div,"r"):
	if "Chr" in l:
		continue
	a=l.split()
	chr,pos,s1a,s2a=a[:4]
	ID=chr+"_"+pos
	## store simulans specific allele
	divhash[ID]=s2a

############################### read input file #########################################################################	
count=0
divhash2=collections.defaultdict(lambda:[])
for sync in filehandle:
	count+=1
	if count%1000000==0:
		print count,"positions processed"
	chr=sync.chr
	pos=sync.pos
	ID=chr+"_"+str(pos)
	pop=sync.subpopulations(pops)
	afl=[]
	if ID in divhash:
		af=0.0
		## go through all populations
		for i in range(len(pop)):
			## sum up allelefrequencies for polymorphic SNPs
			for j in divhash[ID].split("/"):
				if pop[i].cov==0:
					af="NA"
				else:
					if pop[i].counth[j]>=int(options.mincount.split(",")[i]):	
						af+=pop[i].freqh[j]
					else:
						af+=0.0
			afl.append(af)
			divhash2[pops[i]].append(af)
			af=0.0
		
		out.write(chr+"\t"+str(pos)+"\t"+divhash[ID]+"\t"+"\t".join(map(str,afl))+"\n")
		