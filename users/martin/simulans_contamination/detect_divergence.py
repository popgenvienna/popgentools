import modules.RCMH
from modules.CMH import SyncReader,CMHReader
from optparse import OptionParser, OptionGroup
import collections
import math

#########################################################   HELP   #########################################################################
usage="python %prog --input input.sync --output output_cov100.div --species1 1,2,3 --species2 4,5,6 --species1-mincount 1,1,1 --species2-mincount 1,1,1 --min-coverage 50 --max-coverage 200"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
    			D E S C R I P T I O N

This script detects divergent SNPs between two species. Therefore, you may provide multiple populations for each species (--species1 and --species2) which all need to be present in the input sync file (--input). Additionally, you have to provide SNP calling parameters (--min-coverage and --max-coverage), where minimum count has to be set for each species separately (--species1-mincount1 and --species1-mincount2). You can also define the chromosomal arms, which should be used (--chr ; per default: 2L,2LHet,3L,3LHet,3R,3RHet,4,X,Xhet,Yhet). Per default the script is expecting a sync file as an input. However, it also handles CMH outputs (with a P-value at the end). Then, you need to put --type CMH
The output consists of the chromosome, the positions, the allele(s) in species1, the allelel(s) in species2 followed by the synced populations defined with --species1 and --species2

""") 

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--output", dest="output", help="A output file")
parser.add_option("--species1", dest="s1", help="populations of species1 in the format 1,3,4 (according to the position in the sync file)")
parser.add_option("--species2", dest="s2", help="populations of species2 in the format 5,6 (according to the position in the sync file)")
parser.add_option("--species1-mincount",dest="s1m",help="minimum count for the populations of species1 in the format 2,3,3 ")
parser.add_option("--species2-mincount",dest="s2m",help="minimum count for the populations of species1 in the format 2,3 ")
parser.add_option("--min-coverage",dest="mincoverage",help="The minimum coverage")
parser.add_option("--max-coverage",dest="maxcoverage",help="The maximum coverage, if value is below 1, the top xxx % will be removed, if the value is larger than 1, a fixed threshold of this value will be applied to all populations, e.g. 70 for a constant threshold of 70 and 0.02 for a 2% cutoff")
parser.add_option("--type",dest="type",help="file type: cmh output or sync file",default="sync")
parser.add_option("--chr", dest="chr", help="chromosomal arms to be considered in the format: 2L,2LHet,3L,3LHet,3R,3RHet,4,X,Xhet,Yhet",default="2L,2LHet,3L,3LHet,3R,3RHet,4,X,Xhet,Yhet")
parser.add_option_group(group)
(options, args) = parser.parse_args()

########################################################### code ###################################################

def filtercounts(mincount,pophash):
	H=collections.defaultdict(lambda:0)
	for k,v in pophash.items():
		if v >0 and v >=mincount:
			H[k]=v
	return H

############################################## load parameters and open output files ###################################

s1m=map(int,options.s1m.split(","))
s2m=map(int,options.s2m.split(","))
mincoverage=int(options.mincoverage)
sp1=map(int,options.s1.split(","))
sp2=map(int,options.s2.split(","))
ofw=open(options.output,"w")
apm=s1m+s2m
allpop=sp1+sp2

ofw.write("Chr\tPos\tallele_spec1\tallele_spec2\t"+"\t".join(["spec1_pop: "+str(x) for x in sp1])+"\t"+"\t".join(["spec2_pop: "+str(x) for x in sp2])+"\n")
#out2=open(options.output+"_missing_snps","w")

############################################## test for format of input file ###################################

if options.type=="sync":
	filehandle=SyncReader(options.input)
else:
	filehandle=CMHReader(options.input)

############################### max coverage #########################################################################
count=0
if float(options.maxcoverage)<1:
	cutoff=float(options.maxcoverage)
	covdict=collections.defaultdict(lambda:[])
	covlist=[]
	## for each population make a dictionary entry with a list of all coverages
	for sync in filehandle:
		apops=sync.subpopulations(allpop)
		for i in range(len(apops)):
			covdict[i].append(int(apops[i].cov))
		count+=1
		if count%1000000==0:
			print count,"positions processed"
	## calculate max coverage for each population
	for k in sorted(covdict.keys()):
		##make a dictionary of all unique coverages to later detect the next smaller coverage 
		covclass=collections.defaultdict(lambda:0)
		count=1
		for item in sorted(set(covdict[k])):
			covclass[item]=count
			count+=1
		##make a reversed dictionary with the indices as keys and the coverages as values
		revcovclass=dict(zip(covclass.values(),covclass.keys()))
		## determine the class with the top 2% coverages by sorting all coverages from largest to smallest and then slice the list at the position cutoff*length_list (round down to the lower integer)
		covpos=sorted(covdict[k],reverse=True)[int(math.floor(len(covdict[k])*cutoff))]
		##append the next smaller coverage to a list
		covlist.append(revcovclass[covclass[covpos]-1])
		print "population:",allpop[k],", coverage cutoff of top ",cutoff*100,"%:", revcovclass[covclass[covpos]-1]
	maxcov=covlist
	
if float(options.maxcoverage)>1:
	maxcov=int(options.maxcoverage)

############################### calculate divergence #########################################################################

if options.type=="sync":
	filehandle=SyncReader(options.input)
else:
	filehandle=CMHReader(options.input)

count=0
for sync in filehandle:
	count+=1
	if count%1000000==0:
		print count,"positions processed"
	chr=sync.chr
	pos=sync.pos

	# extract the correct populations
	apops=sync.subpopulations(allpop)
	# test if coverage is valid
	if not modules.RCMH.Utility.coverageValid(apops, mincoverage, maxcov):		
		continue
	
	# exclude deletions
	delcount=modules.RCMH.Utility.getDeletionCount(apops)
	if delcount>=1:
		continue
	
	sp1p=sync.subpopulations(sp1)
	sp2p=sync.subpopulations(sp2)
	test=0
	for i in range(len(sp1p)):
		for j in range(len(sp2p)):
			h1=filtercounts(s1m[i],sp1p[i].counth)
			h2=filtercounts(s2m[j],sp2p[j].counth)
			for item in h1.keys():
				if h2.has_key(item):
					test=1
	if test==1:
		continue
	else:
		ofw.write(chr+"\t"+str(pos)+"\t"+"/".join(h1.keys())+"\t"+"/".join(h2.keys())+"\t"+"\t".join([x.__str__() for x in sp1p])+"\t"+"\t".join([x.__str__() for x in sp2p])+"\n")