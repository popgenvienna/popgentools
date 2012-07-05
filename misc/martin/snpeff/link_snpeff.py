import sys
from optparse import OptionParser, OptionGroup
import collections


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
	
t=open(str(options.snpeff),"r")	
dict1=collections.defaultdict(lambda:[])
for l in t:
	if l.rstrip('\n')=='' or "#" in l:
		continue
		l=l.replace(" ","_")
	a=l.split("\t")
	Chromosome=a.pop(0)
	Position=a.pop(0)
	Reference=a.pop(0)
	Change=a.pop(0)
	Changetype=a.pop(0)
	Homozygous=a.pop(0)
	Quality=a.pop(0)
	Coverage=a.pop(0)
	Warnings=a.pop(0)
	Gene_ID=a.pop(0)
	Gene_name=a.pop(0)
	Bio_type=a.pop(0)
	Trancript_ID=a.pop(0)
	Exon_ID=a.pop(0)
	Exon_Rank=a.pop(0)
	Effect=a.pop(0)
	if "UPSTREAM" in  Effect or "DOWNSTREAM" in Effect:
		bases_away=Effect.split(":")[1][1:].replace(" ","_")
	else:
		bases_away="na"
	if "WITHIN" not in Effect and ":" in Effect:
		Effect=Effect.split(":")[0]
	elif "WITHIN" in Effect:
		Effect.replace(" ","_")
	
	old_AA_new_AA=a.pop(0)
	Old_codon_New_codon=a.pop(0)
	Codon_Num=a.pop(0)
	CDS_size=a.pop(0)
	Custom_interval_ID=a.pop(0)
	effectlist=[Reference,Change,Gene_ID,Gene_name,Effect,bases_away,old_AA_new_AA,Old_codon_New_codon]
	#effectlist[effectlist.index("")]="na"
	dict1[Chromosome+"_"+Position].append("\t".join(effectlist))		
	
	
cand=open(str(options.igv),"r")
finallist=[]
for l in cand: 
	if l.rstrip()!="":
		a=l.split("\t")
		ids=a[0]+"_"+a[2]
		if ids in dict1:
			for i in set(dict1[ids]):
				print l.rstrip()+"\t"+i