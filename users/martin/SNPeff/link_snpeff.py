import sys
from optparse import OptionParser, OptionGroup
import collections

#Author: Martin Kapun
#########################################################   HELP   #########################################################################
usage="python %prog --snpeff SNPs.snpeff --input SNP_set.txt > SNP_set.genes"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P :
_________

Description:

This script links the output of SNPeff for all SNPs (-s) to a given set of SNPs (-i). The script appends information from the SNPeff output at the end of each line. Only a certain amount of columns from the SNPeff output is printed (see example). In some cases, SNPs can have more than one effect (e.g. overlapping genes), these SNPs will be printed multiple times for all SNP effects (e.g. nonsynonymous and intergenic). 

Input:

--input:
Any file, which  contains at least two columns with the Chromosome and the position

--snpeff:
The output of SNPeff 2.0.3 for the whole SNPset

Output: 

The script will append the following columns to each SNP in the input:
c1:	Chromosome (of the SNP_set file)
c2:	Position (of the SNP_set file)
c3:	Reference allele (SNP1)
c4:	Alternative allele (SNP2)
c5:	GeneID (e.g. FBgn0200234): If SNP is located within or close to a gene, else \"nan\" is printed.
c6:	Genename (e.g. fru): If SNP is located within or close to a gene, else \"nan\" is printed.
c7:	Effect: (e.g. Intergeneic, Non-synonymous, Synonymous, etc.), see SNPeff homepage for more information.
c8:	Location_to_gene (bases away from next gene): If the SNP is located within the predefined distance of a gene, the distance will be printed in bp, otherwise \"nan\" is printed.
c9:	AA_change: If the SNP is loacted in a coding region, the reference and alternative amino acid will be printed separated by a "/", otherwise \"nan\" is printed.
c10:	codon_change: If the SNP is loacted in a coding region, the reference and alternative codon will be printed separated by a "/", otherwise \"nan\" is printed.
see example:

Chromosome      Position	SNP1	SNP2	GeneID		Genename	Effect  		Location_to_gene	AA_change	codon_change
3R		20573236	C	T	FBgn0004885	tok		INTRON			nan			nan			nan
X		18160684	A	G	nan		nan		INTERGENIC		nan			nan			nan
2L		16742966	A	G	FBgn0032600	CG17912		UPSTREAM		438_bases		nan			nan
3L		15230162	A	T	FBgn0029114	Tollo		SYNONYMOUS_CODING	nan			A/A			gcA/gcT
""")
	
	
#########################################################   CODE   #########################################################################

parser.add_option("-s", "--snpeff", dest="snpeff", help="SNPeff output")
parser.add_option("-i", "--input", dest="input", help="Text file with set of SNPs")
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
	if l.rstrip('\n')=='' or l[0] == "#":
		continue
	a=l.split("\t")
	Chromosome,Position,Reference,Change,Changetype,Homozygous,Quality,Coverage,Warnings,Gene_ID,Gene_name,Bio_type,Trancript_ID,Exon_ID,Exon_Rank,Effect,old_AA_new_AA,Old_codon_New_codon,Codon_Num,Codon_Degeneracy,CDS_size,Codons_around,AAs_around,Custom_interval_ID=a
	if Gene_ID=="":
		Gene_ID="nan"
	if "Gene_" in Gene_ID:
		continue
	if Gene_name=="":
		Gene_name="nan"
	if "UPSTREAM" in  Effect or "DOWNSTREAM" in Effect:
		bases_away=Effect.split(":")[1][1:].replace(" ","_")
	else:
		bases_away="nan"
	if "WITHIN" not in Effect and ":" in Effect:
		Effect=Effect.split(":")[0]
	elif "WITHIN" in Effect:
		Effect.replace(" ","_")
	if old_AA_new_AA=="":
		old_AA_new_AA="nan"
	if Old_codon_New_codon=="":
		Old_codon_New_codon="nan"
	effectlist=[Reference,Change,Gene_ID,Gene_name,Effect,bases_away,old_AA_new_AA,Old_codon_New_codon]
	dict1[Chromosome+"_"+Position].append("\t".join(effectlist))			
cand=open(str(options.input),"r")
finallist=[]
for l in cand: 
	if l.rstrip()!="" and  l[0]!= "#":
		a=l.split("\t")
		ids=a[0]+"_"+a[1]
		if ids in dict1:
			for i in set(dict1[ids]):
				print l.rstrip()+"\t"+i