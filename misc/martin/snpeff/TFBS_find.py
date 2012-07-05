import sys
import operator
import re
import os
from optparse import OptionParser
import collections

parser = OptionParser()
parser.add_option("-i", "--igv", dest="igv", help="IGV file with 3 fst values and index; Needs to be sorted by index")
parser.add_option("-g", "--gff", dest="gff", help="GTF file just with genes (Not EXONS)")
(options, args) = parser.parse_args()

data1=str(options.igv)
#4	521	528	3	11.294139	4.714833	8.67908	521 	0.29166667 	0.16689149 	0.06666667 	17 	10 	14 	0.606927891721 	0.31546777336 	6.85466158757 	0.918646865663
#4	717	717	1	0.046170	0.996290	0.99629	717 	0.01445715 	0.02894239 	0.00468614 	127 	223 	176 	no_biallelic_SNP	no_biallelic_SNP	no_biallelic_SNP	no_biallelic_SNP
#4	810	836	6	4.824219	1.754157	5.595638	810 	0.10788108 	0.13334126 	0.01255887 	94 	74 	52 	0.434931776818 	0.71541311545 	4.83788879744 	34.7687800888
#4	890	890	1	0.252965	1.203085	1.203085	890 	0.02023028 	0.02310894 	0.00141716 	89 	170 	130 	0.470673414466 	0.510335091123 	0.22159288907 	0.240748136157
#4	896	896	1	0.248378	1.198498	1.198498	896 	0.02019019 	0.0216257 	0.00085359 	79 	136 	115 	0.453646483428 	0.478399586467 	0.213637506101 	0.225571154568
#4	905	987	18	39.771016	3.159621	21.338063	975 	0.42188344 	0.40701141 	0.00102794 	50 	63 	57 	0.47980539532 	0.471144735061 	0.3847787497 	0.370933186677
#

data2=str(options.gff)
#2R	REDfly	regulatory_region	21110141	21111300	.	.	.	ID="Kr_CD1"; Dbxref="Flybase:FBgn0001325", "PMID:2114978", "REDfly:RFRC:00000067.001"; Evidence="reporter construct (in vivo)"; Ontology_term="FBbt:00000095","FBbt:00005304";
#2R	REDfly	regulatory_region	21111115	21111706	.	.	.	ID="Kr_HBg0.6HZ"; Dbxref="Flybase:FBgn0001325", "PMID:2114978", "REDfly:RFRC:00000068.001"; Evidence="reporter construct (in vivo)"; Ontology_term="FBbt:00000095","FBbt:00005649";
#2R	REDfly	regulatory_region	21111574	21113281	.	.	.	ID="Kr_NcS1.7HZ"; Dbxref="Flybase:FBgn0001325", "PMID:2114978", "REDfly:RFRC:00000069.001"; Evidence="reporter construct (in vivo)"; Ontology_term="FBbt:00000095","FBbt:00005304","FBbt:00005649";
#2R

bonobo=collections.defaultdict(lambda:[])
gff=open(data2,"r")
entries=[]
tup=()

for line in gff:
	if "##" not in line:
		#print line.split()[8]
		a=line.split()
		#print a
		#print gtfsplit[8]
		var1=re.search(r"ID=\"(.*)\";.*",a[8])
		var2=re.search(r"Dbxref=\"Flybase:(.*)\"",a[9])
		CRM=var1.group(1)
		TF=var2.group(1)
		#print flybase,genename
		#sys.exit()
		tup=(a[0],a[3],a[4],CRM,TF)
		bonobo[a[0]].append(tup)

hesh={}
#print """Chromosome	Start	End	feature	X1.2	X1.3	X1.4	X2.3	X2.4	X3.4	indexNotRounded	av_FST	sumabdiff	av_FSTbetw	start_gene	end_gene	flybase_name	genename	orientation	bp_from_start	bp_from_end"""

igv=open(data1,"r")
for line in igv:
	if "Chromosome" in line: 
		print line.rstrip()+"\tstart\tend\tCRM\tTF\tdistance_from_start\tdistance_from_end\tlength_of_motif"
	else:
		#counter+=1
		count=0
		datasplit=line.split()
		#print datasplit
		for chro, start,end,CRM,TF in bonobo[datasplit[0]]:
			if int(datasplit[2])>=int(start) and int(datasplit[2])<=int(end):
				print line.rstrip()+"\t"+ start+"\t"+ end+"\t"+CRM+"\t"+TF+"\t"+str(int(datasplit[1])-int(start))+"\t"+str(int(datasplit[1])-int(end))+"\t"+str(int(end)-int(start))
#				count+=1
#		if count==0:
#			print line.rstrip()+"\tna\tna\tna\tna\tna\tna\tna"
#			d=datasplit[0]+"-"+datasplit[1]
#			hesh[d]=1
#for line in igv:
#	datasplit=line.split()
#	a=datasplit[0]+"-"+datasplit[1]	
#	if a not in hesh:
#		print line.rstrip()+"na"+"na"+"no_closeby_gene"
#
##