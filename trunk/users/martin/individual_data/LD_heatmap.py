import sys 
from rpy2.robjects import r
import rpy2.robjects as robjects
from optparse import OptionParser, OptionGroup
import collections
import random

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="\npython %prog --input output_file.consensus --ind 2,3,6,10 --subsample 500 --chromosome 2L --output output_2L"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

The purpose of this script is to extract a subsample (e.g. --subsample 0.5) or a defined number of polymorphic SNPs (e.g. --subsample 500) in the Individuals defined with --pop from a consensus file produced with extract_consensus.py (--input). For this SNP-dataset the script will calculate all pairwise r^2 values (according to Hill and Robertson 1968). There will be two output files: The first is a tab delimited file containing a tabular representation of the distance matrix, i.e. the columns consist of Chromosome, Position of SNP1, position of SNP2, r^2, alleles of SNP1 and alleles of SNP2.
The second file is a visual representation of the distance matrix using the LDheatmap R package, which needs to be installed before. This figure will show the physical genomic positions of SNPs used r^2 values highlighted in colors and for chromosomes with cosmopolitan inversion the breakpoints and the name of the inversion. Due to memory constraints, this script can only process one chromosmome at a time, which needs to be defined with (--chromosome).
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="input", help="Input consensus file")
parser.add_option("--output", dest="o", help="output file")
parser.add_option("--subsample", dest="r", help="subsample certain percentage, e.g 0.1, or a defined number of SNPs, e.g. 500",default=1.0)
parser.add_option("--ind", dest="i", help="indivduals used for the analysis")
parser.add_option("--chromosome", dest="c", help="chromosome used for the analysis")

parser.add_option_group(group)
(options, args) = parser.parse_args()


def inversionrect(chromo,steps):
	''' draw lines for inversionbreakpoints in the Heatmap!'''
	incoo=[]
	if chromo=="2L":
		incoo.append((2204115*steps,13212655*steps,"In(2L)t"))
	elif chromo=="3R":
		incoo.append((15922589*steps,28234649*steps,"In(3R)C"))
		incoo.append((17058538*steps,24780914*steps,"In(3R)Mo"))
		incoo.append((12401910*steps,20580518*steps,"In(3R)Payne"))
	elif chromo=="2R":
		incoo.append((11319930*steps,16163246*steps,"In(2R)Ns"))
	elif chromo=="3L":
		incoo.append((3160436*steps,16362055*steps,"In(3L)P"))
	else:
		return "NA"
	return incoo

def doubleallelecodes(line,line1,individuals):
	'''calcualte allele frequencies for each SNP and haplotype frequencies for two SNPs'''
	## get data of first SNP
	chr,pos,ref,mel36,pops,qual,alleles=line.split()
	## get data from second SNP
	chr1,pos1,ref1,mel361,pops1,qual1,alleles1=line1.split()
	codehash=[]
	## loop over inversions and st
	a1,a2,fa=0.0,0.0,0.0
	fu1,fu2,fua="","",[]
	for ind in individuals:
		## get alleles of all indivduals , e.g AAGA
		fu1+=pops[ind]
		fu2+=pops1[ind]
		## get haplotypes of alle individuals, e.g. [AG,AG,AA,...]
		fua.append(pops[ind]+pops1[ind])	
	if "N" not in fu1 and "N" not in fu2:
		## calculate allele frequencies and haplotype frequencies
		a1=fu1.count(fu1[0])/float(len(fu1))
		a2=fu2.count(fu2[0])/float(len(fu2))
		fa=fua.count(fua[0])/float(len(fua))
		codehash=[a1,a2,fa,fu1,fu2]
		#print chr,pos,pos1,inversion,fu1,fu2,fua,a1,a2,fa
	else:
		codehash="NA"
	return codehash

def rsquared(p1,q1,x11):
	''' Hill Robertson 1968'''
	p2=1-p1
	q2=1-q1
	import math
	## calculate D
	D=x11-p1*q1
	if p1*p2*q1*q2==0:
		return "NA"
	else:
		r2=(D**2)/(p1*p2*q1*q2)
		#print p1,p2,q1,q2,x11,r2
		return r2

def execute(snplist,chr):
	out=open(options.o+".dist","w")
	count1,count2,count=0,0,1
	rhash,alist,blist=[],[],[]
	## loop through first SNP
	for i in snplist:
		if count%(float(options.r)/10)==0:
			print count,"positions processed, now at position:",i.rstrip()
		count+=1
		pos=int(i.split()[1])
		## loop through second SNP
		for j in snplist[count1:]:
			pos1=int(j.split()[1])
			## get allele frequencies
			if pos1<=pos:
				continue
			afhash=doubleallelecodes(i,j,individual)
			if afhash=="NA":
				continue
			## calculate r2 from Allele frequencies
			fi,fj,ftotal,allele1,allele2=afhash
			rsq=rsquared(fi,fj,ftotal)
			if rsq=="NA":
				continue
			rhash.append(rsq)
			##store positions of SNP1
			alist.append(pos)
			##store positions of SNP2
			blist.append(pos1)
			#print inversion,rsq,pos,pos1
			out.write("\t".join(map(str,[chr,pos,pos1,rsq,allele1,allele2]))+"\n")
		#print count1,i.rstrip()
		count1+=1
		
	#### make x-axis labels based on the assumption that the labels go from 0-1 and now you need to scale the whole genome relative to these borders
	binlist,labellist=[],[]
	##last SNP pos
	upper=max(alist+blist)
	## first SNP pos
	lower=min(alist+blist)
	## make sure that there is at least ONE SNP
	if upper-lower!=0:
		## here caluclate the relative step of one basepair
		step=1/float(upper-lower)
		invcoo=inversionrect(chr,step)
		## set step size to 2mb
		co=2000000
		## this is the stepsize between the ticks
		stepsize=co*step
		## this is the start position
		start=(co-lower)*step
		## bin the steps in a list
		binlist.append(start)
		labellist.append(str(co/1000000)+"mb")
		co+=2000000
		start+=stepsize
		## append ticks until the step is larger than the position of the last SNP 
		while(co<upper):
			labellist.append(str(co/1000000)+"mb")
			binlist.append(start)
			co+=2000000
			start+=stepsize

		## convert python to R
		cp=robjects.vectors.FloatVector(rhash)
		al=robjects.vectors.IntVector(alist)
		bl=robjects.vectors.IntVector(blist)
		bins=robjects.vectors.FloatVector(binlist)
		labels=robjects.vectors.StrVector(labellist)
		r.assign('values',cp)
		r.assign('al',al)
		r.assign('bl',bl)
		r.assign('bins',bins)
		r.assign('labels',labels)
		## open graphics device and load libraries
		r('library("LDheatmap")')
		r('png("'+options.o+"_"+chr+'.png",width=5000,height=5000)')
		##convert distance list to distance matrix
		r('x.names <- sort(unique(c(al, bl)))')
		r('x.dist <- matrix(0, length(x.names), length(x.names))')
		r('dimnames(x.dist) <- list(x.names, x.names)')
		r('x.ind <- rbind(cbind(match(al, x.names), match(bl, x.names)),cbind(match(bl, x.names), match(al, x.names)))')
		r('x.dist[x.ind] <- rep(values, 2)')
		#print r('t(arev(x.dist))')
		## make LDHeatmap grid object based on the r2 values. Use the topo color palette and put the Chromosome and Inversion in the title. Also print the number of SNPs used.Rotate the whole heatmap by 270 degrees.
		r('M<-LDheatmap(x.dist,sort(unique(c(al, bl)),decreasing=F),color=topo.colors(20),flip=T,geneMapLabelX=10000,title="")')
		## add an X-Axis above heatmap and use the labels generated above
		r('la<-LDheatmap.addGrob(M, grid.xaxis(label=labels,at=bins,main=F,gp=gpar(cex=10),name="axis"),height=0)')
		## add inversion breakpoints
		if invcoo!="NA":
			invcount=0
			alphabet=["a","b","c","d","e","f","g","h"]
			for coord in invcoo:
				#print coord
				## add red line for the inversion boundaries
				r('l'+alphabet[invcount+1]+'<-LDheatmap.addGrob(l'+alphabet[invcount]+', grid.lines(x=c('+str(coord[0])+','+str(coord[1])+'),y='+str(1.1+(invcount/float(5)))+',gp=gpar( lwd=8,col="red")),height='+str(0.1+(invcount/float(500)))+')')
				## add label for the inversion 
				r('l'+alphabet[invcount+2]+'<-LDheatmap.addGrob(l'+alphabet[invcount+1]+', grid.text("'+str(coord[2])+'",x='+str(coord[0])+',y='+str(1.3+(invcount/float(5)))+',gp = gpar(cex = 5)),height='+str(0.1+(invcount/float(500)))+')')
				invcount+=2
		## make everything white.
		r('grid.edit("axis", gp = gpar(col = "white"))')
		## and then just make the ticks and the labels black
		r('grid.edit(gPath("axis", "labels"), gp = gpar(col = "black"))')
		r('grid.edit(gPath("axis", "ticks"), gp = gpar(col = "black",lwd=4))')
		## resize the linewidth of the segments
		r('grid.edit(gPath("geneMap", "segments"), gp = gpar(lwd = 0.2))')
		## increae the size of the color key labels
		r('grid.edit("Key",  gp = gpar(cex = 8))')
		## increase the size of the title
		#r('grid.edit(gPath("heatMap", "title"), gp = gpar(cex=0))')
		
		r('dev.off()')


individual=map(int,options.i.split(","))
chromosome=options.c
#idhash=dict(zip(individual,positions))

fullsnplist=[]
count=1
### read all "proper SNPs" into the memory
for l in open(options.input,"r"):
	if count%10000000==0:
		print count,"positions processed"
	count+=1
	a=l.split()
	if a[0]!=chromosome:
		continue
	## only use positions with more than 1 allele
	if len(a[-1].split("/"))!=2 or "/" in a[3]:
		continue
	## test if no N and polymorphuc
	al=""
	for i in individual:
		al+=a[4][i]
	if "N" in al:
		continue
	if len(set(al))==1:
		continue
	## use proportional subset for calculations. set to 1 to use all
	if float(options.r)<1:
		if random.random() <= float(options.r):
			fullsnplist.append(l)
	## use a defined number of SNPs.
	else:
		fullsnplist.append(l)
if float(options.r)>=1:		
	newdict=random.sample(range(len(fullsnplist)),int(options.r))
	newfullsnps=[]
	for i in sorted(newdict):
		newfullsnps.append(fullsnplist[i])

print "done"
## for each Chromosome loop through all possible combinations of SNPs
if float(options.r)<1:
	execute(fullsnplist,chromosome)
else:
	execute(newfullsnps,chromosome)
	
	

