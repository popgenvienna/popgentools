import sys
from optparse import OptionParser,OptionGroup
from modules.SNPeff import SNPeffReader
import collections
import re


# Authors:
# Martin Kapun
# Robert Kofler

#########################################################   HELP   #########################################################################
print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python snp_summary.py -i SNPs -j candidates -o out -v yes -g yes -f yes 
2)	This script compares two snpeff outputs and produces several summaries. The script requires files with the exxtensions \".genes\" and \".igv\" for each dataset to be located in the same folder with the same name. Therefore, for parameters -i and -j put the candidate files without extensions, the program should find them . The   In the general output (out.summary), it produces tables showing the SNP counts and frequencies for (I) SNP effects and (II) distribution on the chromosomes for both inputs (-i and -j) and (III) the proportion of SNPs on the 3R chromosome located inside or outside the Payne inversion and (IV) the ratio of SNPs on the X and and the autosomes. Optionally it produces tables of genes in file 1 (-g) or file 2 (-f). Additionally the script produces a list of genes which are overlapping between the two input files (-v) and another output file containing the SNPs overlapping in both datasets.
	""")

#########################################################   CODE   #########################################################################

parser.add_option("--snpeff", dest="snpeff", help="input file 1: Needs to be an output from the candidates_link_snpeff.py program")
parser.add_option("--candidates", dest="candidates", help="input file 2: Needs to be an output from the candidates_link_snpeff.py program")
parser.add_option("-o", "--output", dest="p", help="define general output file")
parser.add_option_group(group)
(options, args) = parser.parse_args()

def load_candidates(file):
	canddict=collections.defaultdict(lambda:{})
	for l in open(file,"r"):
		a=re.split("\s+",l)
		chrom=a[0]
		pos=int(a[1])
		canddict[chrom][pos]=1
	return canddict



candhash=load_candidates(options.candidates)



## fill in information in the different statistics
annstats=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
chrstats=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
inverstat=collections.defaultdict(lambda:0)
allcount=0
candscount=0
candgenes=set({})
allgenes=set({})
# parse line and create statistics
for se in SNPeffReader(options.snpeff):
	effect=se.effect
	chrom=se.chr
	pos=se.pos
	annstats[effect]["all"]+=1
	chrstats[chrom]["all"]+=1
	allcount+=1
	allgenes.add((se.geneid,se.genename))
	if((chrom in candhash) and (pos in candhash[chrom])):
		annstats[effect]["cand"]+=1
		chrstats[chrom]["cand"]+=1
		candscount+=1
		candgenes.add((se.geneid,se.genename))
		
	# INVERSION Payne: 3R - 12253902 - 20565305
	if(chrom=="3R" and pos>=12253902 and pos <=20565305):
		inverstat["all"]+=1
		if(chrom in candhash and pos in candhash[chrom]):
			inverstat["cand"]+=1

###
### PRINT
###
	
print "feature\tall SNPs {0}; candidates {1}".format(allcount,candscount)	
# print effect stats
for effect in annstats:
	all=annstats[effect]["all"]
	cand=annstats[effect]["cand"]
	print "feature\t{0}\t{1}\t{2}".format(effect,all,cand)

#chromosome
for chrom in chrstats:
	all=chrstats[chrom]["all"]
	cand=chrstats[chrom]["cand"]
	print "chr\t{0}\t{1}\t{2}".format(chrom,all,cand)
	
#inversion
print "inverstat\t{0}\t{1}".format(inverstat["all"],inverstat["cand"])

for geneid, genename in candgenes:
	print "geneid_cand\t{0}\t{1}".format(geneid,genename)

for geneid, genename in allgenes:
	print "geneid_all\t{0}\t{1}".format(geneid,genename)



