import sys
import collections
from optparse import OptionParser,OptionGroup
import copy
from modules.LDUtil import LDIO,Util
import numpy
from modules.CMH import PopIO

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python LD_dist.py -i candidates.cmh -s full.cmh -b 10 -m 1000 -o out.ld
2)	This script needs the rpy2 package, which can be downloaded from here: http://sourceforge.net/projects/rpy/files/rpy2/
3)	This script takes the set of SNPs in the proximity of candidates SNPs and spits out the average p-values for SNPs in bins until a maximum distance. E.g. -m 500 and -b 20 will result in bins of 50bp length from -500 to +500bp around the candidate.
	""") 
#########################################################   CODE   #########################################################################

parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")
parser.add_option("-s", "--snps", dest="snps", help="cmh output with all SNPs")
parser.add_option("-b", "--bins", dest="bins", help="number of bins")
parser.add_option("-m", "--maxdist", dest="maxdist", help="maximum distance to candidate.")
parser.add_option("-o", "--out", dest="out", help="outputfile for boxplot")
parser.add_option("--measure", dest="measure", help="What should be calculated, median (median), geometric mean (gm)")
parser.add_option_group(group)
(options, args) = parser.parse_args()

# 1: Lade die candidaten SNPs
chrh,candl = LDIO.read_candidatehash(options.inp, int(options.maxdist) )

# 2: Iteriere uber alle SNPs
for l in open(str(options.snps),"r"):
	snp=PopIO.parse_cmhline(l)
	candidates=chrh[snp.chr][snp.pos]
	for cand in candidates:
		cand.appendSNP(snp)

# 3: itererier ueber alle canidaten
ofh=open(options.out,"w")
for cand in candl:
	bins=cand.distributeToBins(int(options.bins))
	toprint=[]
	toprint.append(cand.chr)
	toprint.append(cand.pos)
	toprint.append(cand.pvalue)
	for t in bins:
		value=0
		pvalar=[snp.pvalue for snp in t]
		if options.measure == "median":
			value=Util.median(pvalar)
		elif options.measure=="gm":
			value=Util.geometricmean(pvalar)
		else:
			raise Exception("Unknown option for --measure: "+options.measure)
	
		if value is None:
			value="na"
		toprint.append(value)
	
	s="\t".join(map(str,toprint))
	ofh.write(s+"\n")






#out=open(str(options.out),"w")
#binlist,binhash,bh=[],{},{}
#cand=posh(str(options.inp))
#full=posh(str(options.snps))
#bin=int(options.bins)
#maxdist=int(options.maxdist)
#min1=-maxdist
#for i in range(1,int(bin)+1):
#	max1=min1+(maxdist/int(bin)*2)
#	binlist.append(str(min1))
#	binhash[min1]=[]
#	min1=max1
#binlist.append("cand")	
#binhash["cand"]=[]
#counter=0
#out.write("Chrom\tPos\t"+"\t".join(binlist)+"\t"+"\t".join(binlist)+"\n")
#for chrom,poshash in cand.items():
#	for pos,pval in poshash.items():
#		counter+=1
#		pav,psd,plist=[],[],[]
#		minc,min1=int(pos)-maxdist,-maxdist
#		for i in range(0,int(bin)):
#			maxc=minc+(maxdist/int(bin)*2)
#			max1=min1+(maxdist/int(bin)*2)
#			for snppos,snppval in full[chrom].items():
#				if int(snppos)>=minc and int(snppos)<maxc and int(snppos)!=int(pos):
#					plist.append(-numpy.log10(float(snppval)))
#					binhash[min1].append(-numpy.log10(float(snppval)))
#			if len(plist)>1:
#				pav.append(str(lmean(plist)))
#				psd.append(str(lstdev(plist)))
#			elif len(plist)==1:
#				pav.append(str(plist[0]))
#				psd.append("na")
#			else:
#				pav.append("na")
#				psd.append("na")	
#			plist=[]
#			min1=max1
#			minc=maxc
#		pav.append(str(-numpy.log10(float(pval))))
#		psd.append("0")
#		binhash["cand"].append(-numpy.log10(float(pval)))
#		out.write(str(chrom)+"\t"+str(pos)+"\t"+"\t".join(pav)+"\t"+"\t".join(psd)+"\n")
#		if counter%2==0:
#			print str(counter)+" SNPs processed"
#out2=open(str(options.out)+"_total_Average","w")
#out2.write("bins\tav_pvalues\tSD_pvalues\tcount\n")
#for bin, valuelist in sorted(binhash.items()):
#	out2.write(str(bin)+"\t"+str(lmean(valuelist))+"\t"+str(lstdev(valuelist))+"\t"+str(len(valuelist))+"\n")
