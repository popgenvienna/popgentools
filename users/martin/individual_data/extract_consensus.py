import modules.RCMH
import math
from modules.CMH import SyncReader,CMHReader
import collections
from optparse import OptionParser, OptionGroup
from rpy2.robjects import r
import rpy2.robjects as robjects


#Author: Martin Kapun


#########################################################   HELP   #########################################################################
usage="\npython %prog --input individuals.sync --min-coverage 20 --max-coverage 0.05 --min-count 20 --CI-mode individual --output output_file\n\nor:\npython %prog --input individuals.sync --min-coverage 20 --max-coverage 100,100,100,100,100,90+200,200,200,199,100 --min-count 20 --CI-mode individual --output output_file"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

The purpose of this script is to extract the sire allele from F1 hybrids sequenced as indivduals. The input has to be a sync file with the sequences of the dam (mel36) in the first column, followed by all indivduals. Note that only the Autosomes X,2L,2R,3L,3R and 4 will be considered. In contrast to the individuals, Mel36 has been sequenced as a pool of 10 females, therefore the mean allele frequency is not necessarily 50%. Similarily to the scripts in PopPoolation this script uses cummulative counts of alleles across all synced indivduals/populations to test for minimum count (--min-count). Additionally you have to set a minimum coverage (--min-coverage) threshold, which will be applied to all populations/individuals. 
You have two options for the maximum coverage thresholds (--max-coverage). You can either set the threshold in percent, e.g. 0.05. Then the script will calculate the cutoff based on the top 5 percentile. In contrast to PoPoolation this script calculates these thresholds for each chromosome separately. This is necessary, as male indivduals should only have half the coverage on the X. Then, the script will create an additional output with the cutoff threshold as defined by the input, the IDs ad the actual coverage cutoffs. 
The latter information can be alternatively used as a direct input for the maximum coverage threshold to avoid re-calculation of the percentiles. (In the format: Ind1_2L,
At every position, the script 1) determines, whether the reference (mel36) is within the coverage thresholds and 2) the reference is polymoprhic. I.e. the minor allele is above the min-count threshold and the minor allele frequency (AF) is larger than the lower boundary of the 90% binomial confidence interval (BCI) for a 5% allele frequency at the given coverage.
If the refernce does not fullfill the coverage criteria, the allele from the reference genome is used. If the reference is polyorphic, all indivduals will be set to "N", as the sire allele cannot inferred in this case. 
If the reference is not polymorphic, the individuals are tested for alleles different from the reference: Here you have two choices: 1) to calculate 90% BCI for an allele frequency of 50% for each indivdual separately and test whether the allele frequency is within the range (--CI-mode individual) or 2) to pool all individuals carrying the alternative allele and calculate the 95% BCI for an allele frequency of 50% for the pool and use the alternative allele for every individual if the cummulative AF is within the BCI boundaries (--CI-mode all). Tests have shown that method produces false positives. I would therefore strongly suggest method 1 (which is the default).
The output looks a bit like a pileup file. I have introduced a quality string, which allows to reconstruct, why at certain indivduals an "N" was used rather than an allele. See a description of the strings below:

1: sire allele similar to reference
2: sire allele different from reference
c: coverage too low or too high
r: reference is polymorphic
/: alternative AF outside boundaries of BCI for the particular individual
a: more than two alleles

See a template output below:

2L      6921    G       G       NGGNNNGTGNG     c11ccc121c1     T/G
2L      6922    A       A       AAANNNAAANA     111ccc111c1     A
2L      16933   T       T       TTTNTNTTTNT     111/1c111c1     A/T
2L      19795   T       A/T     NNNNNNNNNNN     rrrrrrrrrrr     A/T
2L      20008   A       A       AGGNANGAANA     122c1c211c1     A/G


Col1: Chromosome
Col2: position
Col3: allele from the reference
Col4: allele(s) from sequenced mel36
Col5: Alleles for each individual; in the same order as in the sync file
Col6: Quality strings for each individual; in the same order as in the sync file
Col7: all alleles above the min-count threshold

In addition, two more outputs will be produced:
1)The file with the extension *.af contains the allele frequencies of individuals with non-reference alleles full-filling all criteria. After columns for chromsome and position, all individuals are listed in the same order as in the sync file. Missing values are indicated as "NA".
2) The file with the extension *_ref.af lists the polymorphic sites in the reference. It contains columns for Chromosome, position, reference genome allele, alleles detected and the frequency of the minor allele.

""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="A synchronized input file")
parser.add_option("--min-coverage", dest="c", help="minimum coverage threshold")
parser.add_option("--min-count", dest="m", help="minimum count threshold")
parser.add_option("--output", dest="o", help="output file(s)")
parser.add_option("--max-coverage", dest="a", help="maximum coverage cutoff")
parser.add_option("--CI-mode", dest="ci", help="choose between 'all' and 'individual', see Help for details, default = 'individual' ",default="individual")
parser.add_option_group(group)
(options, args) = parser.parse_args()


def CI(x,limit):
	''' calulate exact binomial confidence interval'''
	if 0.0 in x:
		return 0.0,0.0
	else:
		cilim=str(float(1-limit)/2)
		n=sum(x)
		X=x[0]
		p=X/float(n)
		df1l=2*(n-X+1)
		df2l=2*X

		Fl=float(r('qf('+cilim+','+str(df1l)+","+str(df2l)+',lower.tail=F)')[0])
		df1u=df2l+2
		df2u=df1l-2

		Fu=float(r('qf('+cilim+','+str(df1u)+","+str(df2u)+',lower.tail=F)')[0])
		### lower limit
		Ll=X/(X+(n-X+1)*Fl)
		### upper limit
		Lu=(X+1)*Fu/(n-X+(X+1)*Fu)
		return Ll,Lu

def run():
	for pop,alleles in sorted(popalleles.items()):
		### test if coverage of individual is above threshold, else append "N"
		if sum(alleles.values())<int(options.c):
			hapallele.append("N")
			qual.append("c")
		else:
			hapallele.append(refal)
			qual.append("1")

def CI_ind():
	for pop,alleles in sorted(popalleles.items()):
		#print pop
		### test if coverage of individual is above threshold, else append "N"
		if sum(alleles.values())<int(options.c) or sum(alleles.values())> covdict[pop][chrom]:
			#print pop,"c"
			hapallele.append("N")
			qual.append("c")
			freqal.append("NA")
		else:
			## test if alleles found in individual is fullfilling the minimum count criteria else delete the allele
			for k in alleles.keys():
				if k not in fullalleles:
					del alleles[k]
			if len(alleles)==2:
				total=sum(alleles.values())
				## test if CI fits the data:
				if CI([total/2,total/2],0.9)[0]<=min(alleles.values())/float(sum(alleles.values())):
					for k in alleles.keys():
						if k!=refal[0]: 
							#print pop,"2"
							hapallele.append(k)
							qual.append("2")
							freqal.append(alleles[k]/float(sum(alleles.values())))
				else:
					#print pop,"2"
					hapallele.append("N")
					qual.append("/") 
					freqal.append("NA")
			elif len(alleles)==1:
					#print pop,"1"
					hapallele.append(alleles.keys()[0])
					qual.append("1")
					freqal.append("NA")
			else:
				#print pop,"a"
				hapallele.append("N")
				qual.append("a")
				freqal.append("NA")

def CI_all():
	### extract the counts of all alleles from populations carrying the minor allele
	testci=collections.defaultdict(lambda:0)
	for pop,alleles in sorted(popalleles.items()):
		if minor in alleles:
			for allele,counts in alleles.items():
				testci[allele]+=counts
	
	if len(testci)==0:
		run()
		
	else:
		fiftyfifty=sum(testci.values())/2
		## test whether the minor allele freq cummulated over all indivduals carrying it is higher than the lower boundary of the 95% binomial CI based on all reads of these individuals and an expected freq of 50%	
		if CI([fiftyfifty,fiftyfifty],0.95)[0]<=min(testci.values())/float(sum(testci.values())):
		
			for pop,alleles in sorted(popalleles.items()):
				### test if coverage of individual is above threshold, else append "N"
				if sum(alleles.values())<=int(options.c) or sum(alleles.values())>= covdict[pop][chrom]:
					hapallele.append("N")
					qual.append("c")
					freqal.append("NA")
				else:
					for k in alleles.keys():
						if k not in fullalleles:
							del alleles[k]
					if len(alleles)==2:
						for k in alleles.keys():
							if k!=refal[0]:
								hapallele.append(k)
								freqal.append(alleles[k]/float(sum(alleles.values())))
								qual.append("2")
					elif len(alleles)==1:
							hapallele.append(alleles.keys()[0])
							qual.append("1")
							freqal.append("NA")
					else:
						hapallele.append("N")
						qual.append("a")
						freqal.append("NA")
		else:
			run()

############################### max coverage #################################################

#individual=["mel36",52,53,80,100,117,132,136,150,168,85,106]
#iupac={"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["C","T"]}
synccode=["A","T","C","G"]
chromo=["2L","2R","3L","3R","4","X"]
popstotest=range(3,len(open(options.input).readline().split()[3:])+3)

if "," not in options.a:
	out_cov=open(options.o+".cov","w")
	cutoff=float(options.a)
	covdict_full=collections.defaultdict(lambda:collections.defaultdict(lambda:[]))
	covdict=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
	count=0
	for line in open(options.input,"r"):
		a=line.split()
		chrom=a[0]
		for pop in popstotest:
			if a[pop]!="-":
				#print line[pop], line, pop
				covdict_full[pop][chrom].append(sum(map(int,a[pop].split(":"))))
	
		count+=1
		if count%1000000==0:
			print count,"positions processed"
	string=[]
	chromlist=[]
	for pop,hash1 in sorted(covdict_full.items()):
		string_chrom=[]
		chrom_ID=[]
		for chrom,values in sorted(hash1.items()):
			##make a dictionary of all unique coverages to later detect the next smaller coverage 
			covclass=collections.defaultdict(lambda:0)
			count=1
			for item in sorted(set(values)):
				covclass[item]=count
				count+=1
			##make a reversed dictionary with the indices as keys and the coverages as values
			revcovclass=dict(zip(covclass.values(),covclass.keys()))
			## determine the class with the top 2% coverages by sorting all coverages from largest to smallest and then slice the list at the position cutoff*length_list (round down to the lower integer)
			covpos=sorted(values,reverse=True)[int(math.floor(len(values)*cutoff))]
	
			##append the next smaller coverage to a list
			if covclass[covpos]<2:
				covdict[pop][chrom]=0
			else:
				covdict[pop][chrom]=revcovclass[covclass[covpos]-1]
			string_chrom.append(covdict[pop][chrom])
			chrom_ID.append(str(pop)+"_"+chrom)
		string.append(string_chrom)
		chromlist.append(chrom_ID)
	## write the ID's and coverages to an output file. the last line can be used as the input for another run if you do not want to calucalte everything again	
	out_cov.write(" coverage cutoff of top ",cutoff*100,"%\n")
	out_cov.write("+".join([",".join(x) for x in chromlist ])+"\n")
	out_cov.write("+".join([",".join(map(str,x)) for x in string ])+"\n")
	
	covdict_full=0

else:
	covdict=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
	for i in range(len(options.a.split("+"))):
		covdict[popstotest[i]]=dict(zip(chromo,map(int,options.a.split("+")[i].split(","))))
#print covdict
############################### determine haplotype #########################################

out=open(options.o+".consensus","w")
out2=open(options.o+".af","w")
out3=open(options.o+"_ref.af","w")
count=1
for line in open(options.input,"r"):
	a=line.split()
	chrom=a[0]
	if count%1000000==0:	
		print count,"lines processed"
	count+=1
	fullalleles=collections.defaultdict(lambda:0)
	popalleles=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
	
	## extract the alleles with a count > 0 cummulatively for all populations (fullalleles) and for each population seperatly (popalleles)
	for pop in popstotest:
		popalleles[pop]
		
		### go through all nucleotides per populations and test whether larger than one, if yes, keep!
		for i in range(len(synccode)):
			if int(a[pop].split(":")[i])>0:
			
				## store counts cummulatively for all populations
				fullalleles[synccode[i]]+=int(a[pop].split(":")[i])
				
				## store counts for each population separately
				popalleles[pop][synccode[i]]+=int(a[pop].split(":")[i])
	
	## if an allele is below the minimum threshold defined in the commandline, remove this allele from the cummulative dictionary
	for k,v in fullalleles.items():
		if v < int(options.m):
			del fullalleles[k]

		
	## if overall coverage is too low just print N's for all pops and continue
	if len(fullalleles)==0:
		out.write("\t".join(a[:3])+"\t"+a[2]+"\t"+"N"*(len(popstotest)-1)+"\t"+"c"*(len(popstotest)-1)+"\t-\n")
		continue
	
	############################# extract reference allele ###################################
	
	reference=popalleles[3]
	
	## delete alleles not fullfilling the overall criteria
	for allele in reference.keys():
		if allele not in fullalleles:
			del reference[allele]
			
	## test if reference is above minimum and below maximum coverage 		
	if sum(reference.values())<int(options.c) or sum(reference.values())>covdict[3][chrom]:
		refal=a[2]
		del popalleles[3]
	
	## print "N's"	and exit, if reference is ambiguous
	elif len(reference)>1:
		total=sum(reference.values())
		
		## test if minor allele is within the range of 90% interval of a minimum AF of 0.05 (because reference is a pool of 10 females). If yes, print "N's" and exit, if reference is ambiguous
		if CI([total*0.05,total*0.95],0.9)[0]<=min(reference.values())/float(total):
			out.write("\t".join(a[:3])+"\t"+"/".join(reference.keys())+"\t"+"N"*(len(popstotest)-1)+"\t"+"r"*(len(popstotest)-1)+"\t"+"/".join(fullalleles.keys())+"\n")
			out3.write("\t".join(a[:3])+"\t"+"/".join(reference.keys())+"\t"+str(min(reference.values())/float(total))+"\n")
			continue
			
		## elses use major allele as reference
		else:
			tar= dict((value,key) for key,value in reference.items())
			del popalleles[3]
			refal=tar[max(tar.keys())]
	## if (for whatever reason) there is no allele, take the old reference 
	elif len(reference)==0:
		refal=a[2]
		del popalleles[3]
	## else use the only allele as the new reference allele 	
	else:
		refal=reference.keys()[0]
		del popalleles[3]
	
	### keep two major alleles if more than three alleles in whole samples
	if len(fullalleles)>=2:
		tar= [(key,value) for key,value in fullalleles.items()]
		tar.sort(key=lambda x:-x[1])
		fullalleles=collections.defaultdict(lambda:0)
		for k,v in tar[:3]:
			fullalleles[k]=v
		
	################# determine the sire allele in the indivduals ###########################
	
	freqal=[]  ## list of allelesfreqs from the haplotypes
	hapallele=[] ## list of alleles from the haplotypes
	qual=[] ## list of description symbols for each individual
	if len(fullalleles)==2:
		minor=tar[1][0]
		if options.ci=="individual":
			CI_ind()
		else: 
			CI_all()
	else:
		run()

	out.write("\t".join(a[:3])+"\t"+"/".join(refal)+"\t"+"".join(hapallele)+"\t"+"".join(qual)+"\t"+"/".join(fullalleles.keys())+"\n")
	if len(freqal)!=0 and list(set(freqal))!=["NA"]:
		out2.write("\t".join(a[:2])+"\t"+"\t".join(map(str,freqal))+"\n")
	
