import sys
from optparse import OptionParser, OptionGroup


#Author: Martin Kapun
#version 1.0

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python link_snpeff.py -s SNPs.snpeff -i candidates.igv > candidates.genes
2)	This script links the output of SNPeff for all SNPs (-s) to a given set of SNPS	Steps (-s). The script appends information from the SNPeff output at the end of each line. Only a certain amount of columns from the SNPeff output is finallist.append(ed. (See example)) As in some cases, SNPs can have more than one effect (e.g. If the SNPs is located in the end of one long transcript but not in the CDS of another transcript), these SNPs will be finallist.append(ed multiple times for all SNP effects (e.g. nonsynonymous and intergenic). 
3) When the *.igv file contains a header (starting with \"Chromosome\") the header will also be finallist.append(ed in the output followed by the header for the SNPeff part. Otherwise the script finallist.append(s \"na\" as header for the columns of the candidate input followed by the header form the SNPeff output. See example: 
 
Chromosome      Start   End     feature X1.2    X1.3    X2.3    indexNotRounded av_FST  sumabdiff       av_FSTbetw      Reference       Change  Gene_ID Gene_name       Effect  bases_away      old_AA/new_AA   Old_codon/New_codon
2L      484     485     fst_single      0.01212599      0.00236296      0.00489431      0.277331591697175       0.007244475     0.00976303      0.00489431      G       A       na      na      INTERGENIC      na      na      na
2L      501     502     fst_single      0.00755038      0.00191656      0.00501386      0.0410897803593715      0.00473347      0.00563382      0.00501386      T       A       na      na      INTERGENIC      na      na      na
2L      502     503     fst_single      0.00286019      0.00594038      0.00630994      -0.088467361969693      0.004400285     0.00308019      0.00630994      T       G       na      na      INTERGENIC      na      na      na


c1:	col. 1 of the candidate file 
.
.
.
.
cX:	col. X of the candidate file

cX1:	Reference allele
cX2:	Alternative allele
cX3:	GeneID (e.g. FBgn0200234), if SNP is located within or close to a gene, else \"na\" is finallist.append(ed
cX4:	Genename (e.g. fru), if SNP is located within or close to a gene, else \"na\" is finallist.append(ed
cX5:	SNP-effect: (e.g. Intergeneic, Non-synonymous, Synonymous, etc.)
cX6:	Bases away from next gene: if the SNP is located within the predefined distance of a gene, the distance will be finallist.append(ed in bp, otherwise \"na\" is finallist.append(ed
CX7:	AA change: If the SNP is loacted in a coding region, the reference and alternative AA will be finallist.append(ed separated by a "/", otherwise \"na\" is finallist.append(ed
CX8:	Cododn Change: 	If the SNP is loacted in a coding region, the reference and alternative CODON will be finallist.append(ed separated by a "/", otherwise \"na\" is finallist.append(ed
	""")
	
	
#########################################################   CODE   #########################################################################

parser.add_option("-s", "--snpeff", dest="snpeff", help="SNPeff output")
parser.add_option("-i", "--igv", dest="igv", help="IGV file with candidate SNPs")
parser.add_option_group(group)
(options, args) = parser.parse_args()

#finallist.append( """Chromosome	Start	End	feature	X1.2	X1.3	X1.4	X2.3	X2.4	X3.4	indexNotRounded	av_FST	sumabdiff	av_FSTbetw	Reference	Change	Gene_ID	Gene_name	Effect	bases_away	old_AA/new_AA	Old_codon/New_codon	"""
t=open(str(options.snpeff),"r")
# Chromo	Position	Reference	Change	Change_type	Homozygous	Quality	Coverage	Warnings	Gene_ID	Gene_name	Bio_type	Trancript_ID	Exon_ID	Exon_Rank	Effect	old_AA/new_AA	Old_codon/New_codon	Codon_Num(CDS)	CDS_size	Custom_interval_ID
#2L	274948	C	T	SNP	Hom				FBgn0086855	CG17078		FBtr0300804	FBgn0086855:1	1	SYNONYMOUS_STOP	*/*	TGA/TAA	611	1833	
#2L	7792	A	C	SNP	Hom				FBgn0031208	CG11023		FBtr0300690	FBgn0031208:1	1	NON_SYNONYMOUS_CODING	D/A	GAC/GCC	38	1443
#2L	103762	A	C	SNP	Hom				FBgn0031217	CG11377		FBtr0078104	FBgn0031217:3	3	STOP_LOST	*/S	TAA/TCA	375	1125	
#2L	277646	T	C	SNP	Hom				FBgn0003444	smo		FBtr0078129	CG11561:1	1	START_GAINED: CTG, 5PRIME_UTR: 180 bases from TSS					
#2L	281613	T	C	SNP	Hom				FBgn0003444	smo		FBtr0078129	CG11561:6	6	3PRIME_UTR: 99 bases from transcript end					
#2L	282165	T	A	SNP	Hom				FBgn0003444	smo		FBtr0078129	CG11561:6	6	SPLICE_SITE					
#2L	282523	C	T	SNP	Hom				FBgn0031244	CG11601		FBtr0078137			INTRON				828	
#2L	478029	A	G	SNP	Hom				FBgn0043364	cbt		FBtr0078084	CG4427:3	3	START_LOST	M/T	ATG/ACG	1	1044	
#2L	479096	T	C	SNP	Hom				FBgn0043364	cbt		FBtr0078084			UPSTREAM: 152 bases					
#2L	479958	A	T	SNP	Hom										INTERGENIC					
#2L	540615	G	T	SNP	Hom				FBgn0003963	ush		FBtr0078063			DOWNSTREAM: 74 bases					
#2L	541774	A	G	SNP	Hom				FBgn0010602	lwr		FBtr0078082	CG3018:3	3	SYNONYMOUS_CODING	T/T	ACT/ACC	158	480
#2L	542268	T	G	SNP	Hom				FBgn0010602	lwr		FBtr0078082	CG3018:3	3	5PRIME_UTR: 21 bases from TSS					
#2L	546470	C	T	SNP	Hom				FBgn0031261	nAcRbeta-21C		FBtr0078065	CG11822:7	7	STOP_GAINED	Q/*	CAG/TAG	322	984	
	
	
dict1={}
dict2={}
if t.readline().rstrip('\n')!='' or "#" not in t.readline():
	b=t.readline().split("\t")[0]+"_"+t.readline().split("\t")[1]
	dict1[t.readline().split("\t")[15]]="@".join(t.readline().split("\t")).rstrip("\n")
for line in t:
	if line.rstrip('\n')=='' or "#" in line:
		continue
	l=line.split("\t")

	if l[0]+"_"+l[1]==b:
		dict1[l[15]]="@".join(l).rstrip("\n")
		b=l[0]+"_"+l[1]

	else: 
		dict2[b]=dict1
		b=l[0]+"_"+l[1]
		dict1={}
		dict1[l[15]]="@".join(l).rstrip("\n")

cand=open(str(options.igv),"r")
a=cand.readline()
if "Chromosome" in a:
	print a.rstrip()+"""	Reference	Change	Gene_ID	Gene_name	Effect	bases_away	old_AA/new_AA	Old_codon/New_codon"""
else:
	print "na\t"*len(a.split())+"""Reference	Change	Gene_ID	Gene_name	Effect	bases_away	old_AA/new_AA	Old_codon/New_codon"""
cand=open(str(options.igv),"r")

finallist=[]
for line in cand: 
	if line.rstrip()!="":
		l=line.split("\t")
		if l[0]+"_"+l[2] in dict2:
			for a in dict2[l[0]+"_"+l[2]].items():
				i=a[1].split("@")
				if "SYNONYMOUS_STOP" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15]+"\tna\t"+i[16]+"\t"+i[17])
				if "STOP_LOST" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15]+"\tna\t"+i[16]+"\t"+i[17])
				if "NON_SYNONYMOUS_CODING" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15]+"\tna\t"+i[16]+"\t"+i[17])
				if "START_GAINED" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15].split(":")[0]+"\t"+i[15].split(":")[1]+": "+i[15].split(":")[2]+"\tna\tna")
				if "3PRIME_UTR" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15].split(":")[0]+"\t"+i[15].split(":")[1]+"\tna\tna")
				if "SPLICE_SITE" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15]+"\tna\tna\tna")
				if "INTRON" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15]+"\tna\tna\tna")
				if "START_LOST" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15]+"\tna\t"+i[16]+"\t"+i[17])
				if "UPSTREAM" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15].split(":")[0]+"\t"+i[15].split(":")[1]+"\tna\tna")
				if "INTERGENIC" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\tna\tna\t"+i[15]+"\tna\tna\tna")
				if "DOWNSTREAM" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15].split(":")[0]+"\t"+i[15].split(":")[1]+"\tna\tna")
				if "SYNONYMOUS_CODING" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15]+"\tna\t"+i[16]+"\t"+i[17])
				if "5PRIME_UTR" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15].split(":")[0]+"\t"+i[15].split(":")[1]+"\tna\tna")
				if "STOP_GAINED" in i[15]:
					finallist.append( line.rstrip()+"\t"+i[2]+"\t"+i[3]+"\t"+i[9]+"\t"+i[10]+"\t"+i[15]+"\tna\t"+i[16]+"\t"+i[17])
for line in sorted(set(finallist)):
	print line
