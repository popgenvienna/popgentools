#ex2gff3.py
#Author: Martin Kapun
#2011-01-11
#version 2.2 (for Viola)

from optparse import OptionParser, OptionGroup

#########################################################   HELP   #########################################################################
parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: %prog --input output.exonerate > output.gff3
2)	translates vulgar string from an exonerate output into gff3 file format. The following features are provided in the output gff3: gene, mRNA, exon, CDS, match_part, Indel, Splice Site, Frameshift, Gap	""")

#########################################################   CODE   #########################################################################

parser.add_option("-i", "--input", dest="input", help="Output from exonerate which contains a vulgar string (i.e. --showvulgar set to TRUE)")
parser.add_option_group(group)
(options, args) = parser.parse_args()

exonerate=open(str(options.input),"r")

def shift(self):
	item = self[0]
	del self[:1]
	return item
def reverse(row):
    """func that reverse a list not in place"""
    row.reverse()
    return row
print """##gff-version 3"""
for l in exonerate:
	tup=()
	gffdict={}
	if "vulgar" in l:
		ls=l.split()
		query=l.split()[1]
		starttg=str(int(l.split()[6])+1)
		endtg=str(int(l.split()[7]))
		target=l.split()[5]
		orient=l.split()[8]
		exon1=l.split()[10]
		nuc1=int(starttg)+int(l.split()[12])
		

		gffset=l.split()[10:]
		count=0
		if orient=="+":
			ff=""
			ef=""
			print target+"\tprotein2genome\tgene\t"+starttg+"\t"+str(endtg)+"\t.\t"+orient+"\t.\tID="+query+";Name="+query
			print target+"\tprotein2genome\tmRNA\t"+starttg+"\t"+str(endtg)+"\t.\t"+orient+"\t.\tID="+query+"-R;Name="+query+"-R;Parent="+query
			nuc1=int(starttg)
			while 0<len(gffset):
				field=shift(gffset)
				aa=shift(gffset)
				nuc=shift(gffset)
				count+=1
				tup=(field,nuc,aa)
				gffdict[str(count)]=tup
			gffrange=len(gffdict)
			#print gffdict
			mcount=1
			icount=1
			for i in range(gffrange): 
				feature=gffdict[str(i+1)][0]
				bplength=gffdict[str(i+1)][1]
				codon=gffdict[str(i+1)][2]
				if feature=="M" and ff=="":
					es=nuc1
					ex=nuc1
					nuc1+=int(bplength)
					ee=nuc1
					readingframe=0
					ff="M"
					ef="M"
				if feature=="M" and ff=="F":
					es=nuc1
					nuc1+=int(bplength)
					ee=nuc1
					ff="M"
				if feature=="M" and ff=="G":
					es=nuc1
					nuc1+=int(bplength)
					ee=nuc1
					ff="M"
				if feature=="S" and ff=="M":
					nuc1+=int(bplength)
					ee=nuc1
					ff="MS"
				if feature=="5" and ff=="MS":
					print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
					print target+"\tprotein2genome\texon\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R"
					print target+"\tprotein2genome\tCDS\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R"
					print target+"\tprotein2genome\tsplice_site\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R"	
					nuc1+=int(bplength)
					ff=""
					ef=""	
					mcount+=1
				if feature=="I":
					print target+"\tprotein2genome\tintron\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_intron:"+str(icount)+";Name="+query+"_intron:"+str(icount)+";Parent="+query+"-R"
					nuc1+=int(bplength)
					ff=""
					icount+=1
				if feature=="3" and ff=="":
					print target+"\tprotein2genome\tsplice_site\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff=""
				if feature=="S" and ff=="":
					es=nuc1
					ex=nuc1
					nuc1+=int(bplength)
					readingframe=int(bplength)
					ff="S"	
					ef="M"
				if feature=="M" and ff=="S":
					nuc1+=int(bplength)
					ee=nuc1
					ff="SM"	
				if feature=="S" and ff=="SM":
					nuc1+=int(bplength)
					ee=nuc1
					ff="SMS"	
				if feature=="5" and ff=="SMS":
					print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
					print target+"\tprotein2genome\texon\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R"
					print target+"\tprotein2genome\tCDS\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R"	
					print target+"\tprotein2genome\tsplice_site\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff=""
					ef=""	
					mcount+=1
				if feature=="5" and ff=="SM":
					print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
					print target+"\tprotein2genome\texon\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R"
					print target+"\tprotein2genome\tCDS\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R"	
					print target+"\tprotein2genome\tsplice_site\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff=""
					ef=""	
					mcount+=1					
				if feature=="5" and ff=="MS":
					print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
					print target+"\tprotein2genome\texon\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R"
					print target+"\tprotein2genome\tCDS\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R"	
					print target+"\tprotein2genome\tsplice_site\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff=""
					ef=""	
					mcount+=1
				if feature=="5" and ff=="M":
					print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
					print target+"\tprotein2genome\texon\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R"
					print target+"\tprotein2genome\tCDS\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R"	
					print target+"\tprotein2genome\tsplice_site\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff=""
					ef=""	
					mcount+=1	
				if feature=="5" and ff=="G":
					print target+"\tprotein2genome\tsplice_site\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff=""
				if feature=="5" and ff=="F":
					print target+"\tprotein2genome\tsplice_site\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff=""	
				if  feature=="G" and ff=="M":
					if int(codon)!=0:
						print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
						print target+"\tprotein2genome\tindel\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Name="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Parent="+query+"-R"
						nuc1+=int(bplength)	
						ff="G"
					if int(codon)==0:
						print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
						print target+"\tprotein2genome\tindel\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(bplength))+"_bp;Name="+query+"_indel_len_"+str(int(bplength))+"_bp;Parent="+query+"-R"
						nuc1+=int(bplength)	
						ff="G"
				if feature=="G" and ff=="SM":
					if int(codon)!=0:
						print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
						print target+"\tprotein2genome\tindel\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Name="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Parent="+query+"-R"
						nuc1+=int(bplength)
						ff="G"
					if int(codon)==0:
						print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
						print target+"\tprotein2genome\tindel\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(bplength))+"_bp;Name="+query+"_indel_len-"+str(int(bplength))+"_bp;Parent="+query+"-R"
						nuc1+=int(bplength)
						ff="G"
				if feature=="G" and ff=="F":
					if int(codon)!=0:
						print target+"\tprotein2genome\tindel\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Name="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Parent="+query+"-R"
						nuc1+=int(bplength)
						ff="G"
					if int(codon)==0:
						print target+"\tprotein2genome\tindel\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(bplength))+"_bp;Name="+query+"_indel_len-"+str(int(bplength))+"_bp;Parent="+query+"-R"
						nuc1+=int(bplength)
						ff="G"
				if feature=="G" and ff=="":
					if int(codon)!=0:
						print target+"\tprotein2genome\tindel\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Name="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Parent="+query+"-R"
						nuc1+=int(bplength)
						ff="G"
					if int(codon)==0:
						print target+"\tprotein2genome\tindel\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(bplength))+"_bp;Name="+query+"_indel_len-"+str(int(bplength))+"_bp;Parent="+query+"-R"
						nuc1+=int(bplength)
						ff="G"
				if feature=="F" and ff=="M":
					print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
					print target+"\tprotein2genome\tframeshift\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_frameshift;Name="+query+"_frameshift;Parent="+query+"-R"
					nuc1+=int(bplength)	
					ff="F"
				if feature=="F" and ff=="SM":
					print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
					print target+"\tprotein2genome\tframeshift\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_frameshift;Name="+query+"_frameshift;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff="F"
				if feature=="F" and ff=="G":
					print target+"\tprotein2genome\tframeshift\t"+str(nuc1)+"\t"+str(nuc1+int(bplength)-1)+"\t.\t"+orient+"\t.\tID="+query+"_frameshift;Name="+query+"_frameshift;Parent="+query+"-R"
					nuc1+=int(bplength)
					ff="F"
				if feature=="N" and ff=="S": 
					nuc1+=int(bplength)
					ff="S"
				if feature=="N" and ff=="M": 
					nuc1+=int(bplength)
					ff="M"
				if feature=="N" and ff=="SM": 
					nuc1+=int(bplength)
					ff="SM"
			print target+"\tprotein2genome\tmatch_part\t"+str(es)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R"
			print target+"\tprotein2genome\texon\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R"
			print target+"\tprotein2genome\tCDS\t"+str(ex)+"\t"+str(ee-1)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R"
		
		if orient=="-":
			invert1=[]
			ff=""
			ef=""
			print target+"\tprotein2genome\tgene\t"+str(int(endtg)+1)+"\t"+str(int(starttg)-1)+"\t.\t"+orient+"\t.\tID="+query+";Name="+query
			print target+"\tprotein2genome\tmRNA\t"+str(int(endtg)+1)+"\t"+str(int(starttg)-1)+"\t.\t"+orient+"\t.\tID="+query+"-R;Name="+query+"-R;Parent="+query
			nuc1=int(starttg)-1
			while 0<len(gffset):
				field=shift(gffset)
				aa=shift(gffset)
				nuc=shift(gffset)
				count+=1
				tup=(field,nuc,aa)
				gffdict[str(count)]=tup
			gffrange=len(gffdict)
			#invert1.append(gffdict
			mcount=1
			icount=1
			for i in range(gffrange): 
				feature=gffdict[str(i+1)][0]
				bplength=gffdict[str(i+1)][1]
				codon=gffdict[str(i+1)][2]
				if feature=="M" and ff=="G":
					es=nuc1
					nuc1-=int(bplength)
					ee=nuc1
					ff="M"
					
				if feature=="M" and ff=="F":
					es=nuc1
					nuc1-=int(bplength)
					ee=nuc1
					ff="M"
					
				if feature=="M" and ff=="":
					es=nuc1
					ex=nuc1
					nuc1-=int(bplength)
					ee=nuc1
					readingframe=0
					ff="M"
					ef="M"
					
				if feature=="S" and ff=="M":
					nuc1-=int(bplength)
					ee=nuc1
					ff="MS"
					
				if feature=="5" and ff=="MS":
					invert1.append(target+"\tprotein2genome\texon\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tCDS\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tsplice_site\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					ef=""
					mcount+=1	
					
				if feature=="I":
					invert1.append(target+"\tprotein2genome\tintron\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_intron:"+str(icount)+";Name="+query+"_intron:"+str(icount)+";Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					icount+=1
					
				if feature=="3" and ff=="":
					invert1.append(target+"\tprotein2genome\tsplice_site\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					
				if feature=="S" and ff=="":
					es=nuc1
					ex=nuc1
					nuc1-=int(bplength)
					readingframe=int(bplength)
					
					ff="S"
					ef="M"	
				if feature=="M" and ff=="S":
					nuc1-=int(bplength)
					ee=nuc1
					ff="SM"	
					
				if feature=="S" and ff=="SM":
					nuc1-=int(bplength)
					ee=nuc1
					ff="SMS"	
					
				if feature=="5" and ff=="SMS":
					invert1.append(target+"\tprotein2genome\texon\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tCDS\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tsplice_site\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					ef=""
					mcount+=1
					
				if feature=="5" and ff=="G":
					invert1.append(target+"\tprotein2genome\tsplice_site\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					
				if feature=="5" and ff=="F":
					invert1.append(target+"\tprotein2genome\tsplice_site\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					
				if feature=="5" and ff=="SM":
					invert1.append(target+"\tprotein2genome\texon\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tCDS\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tsplice_site\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					ef=""
					mcount+=1		
								
				if feature=="5" and ff=="MS":
					invert1.append(target+"\tprotein2genome\texon\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tCDS\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tsplice_site\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					ef=""
					mcount+=1	
					
				if feature=="5" and ff=="M":
					invert1.append(target+"\tprotein2genome\texon\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tCDS\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tsplice_site\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_splice;Name="+query+"_splice_site;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff=""
					ef=""
					mcount+=1	
					
				if  feature=="G" and ff=="M":
					if int(codon)!=0:
						invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
						invert1.append(target+"\tprotein2genome\tindel\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Name="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Parent="+query+"-R")
						nuc1-=int(bplength)	
						ff="G"
						
					if int(codon)==0:
						invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
						invert1.append(target+"\tprotein2genome\tindel\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len-"+str(int(bplength))+"_bp;Name="+query+"_indel_len-"+str(int(bplength))+"_bp;Parent="+query+"-R")
						nuc1-=int(bplength)	
						ff="G"
						
				if feature=="G" and ff=="SM":
					if int(codon)!=0:
						invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
						invert1.append(target+"\tprotein2genome\tindel\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Name="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Parent="+query+"-R")
						nuc1-=int(bplength)
						ff="G"
						
					if int(codon)==0:
						invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
						invert1.append(target+"\tprotein2genome\tindel\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(bplength))+"_bp;Name="+query+"_indel_len_"+str(int(bplength))+"_bp;Parent="+query+"-R")
						nuc1-=int(bplength)
						ff="G"
				if feature=="G" and ff=="":
					if int(codon)!=0:
						invert1.append(target+"\tprotein2genome\tindel\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Name="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Parent="+query+"-R")
						nuc1-=int(bplength)
						ff="G"
						
					if int(codon)==0:
						invert1.append(target+"\tprotein2genome\tindel\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(bplength))+"_bp;Name="+query+"_indel_len_"+str(int(bplength))+"_bp;Parent="+query+"-R")
						nuc1-=int(bplength)
						ff="G"
						
				if feature=="G" and ff=="F":
					if int(codon)!=0:
						invert1.append(target+"\tprotein2genome\tindel\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Name="+query+"_indel_len_"+str(int(codon)*3)+"_bp;Parent="+query+"-R")
						nuc1-=int(bplength)
						ff="G"
						
					if int(codon)==0:
						invert1.append(target+"\tprotein2genome\tindel\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_indel_len_"+str(int(bplength))+"_bp;Name="+query+"_indel_len_"+str(int(bplength))+"_bp;Parent="+query+"-R")
						nuc1-=int(bplength)
						ff="G"
						
				if feature=="F" and ff=="M":
					invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tframeshift\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_frameshift;Name="+query+"_frameshift;Parent="+query+"-R")
					nuc1-=int(bplength)	
					ff="F"
					
				if feature=="F" and ff=="SM":
					invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
					invert1.append(target+"\tprotein2genome\tframeshift\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_frameshift;Name="+query+"_frameshift;Parent="+query+"-R")
					nuc1-=int(bplength)
					ff="F"
					
				if feature=="F" and ff=="":
					invert1.append(target+"\tprotein2genome\tframeshift\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_frameshift;Name="+query+"_frameshift;Parent="+query+"-R")
					nuc1-=int(bplength)	
					ff="F"
					
				if feature=="F" and ff=="G":
					invert1.append(target+"\tprotein2genome\tframeshift\t"+str(nuc1-int(bplength)+1)+"\t"+str(nuc1)+"\t.\t"+orient+"\t.\tID="+query+"_frameshift;Name="+query+"_frameshift;Parent="+query+"-R")
					nuc1-=int(bplength)	
					ff="F"
					
				if feature=="N" and ff=="S": 
					nuc1-=int(bplength)
					ff="S"
					
				if feature=="N" and ff=="M": 
					nuc1-=int(bplength)
					ff="M"
					
				if feature=="N" and ff=="SM": 
					nuc1-=int(bplength)
					ff="SM"
					
			invert1.append(target+"\tprotein2genome\texon\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t.\tID="+query+":"+str(mcount)+";Name="+query+":"+str(mcount)+";Parent="+query+"-R")
			invert1.append(target+"\tprotein2genome\tCDS\t"+str(ee+1)+"\t"+str(ex)+"\t.\t"+orient+"\t"+str(readingframe)+"\tID="+query+"_cds;Name="+query+"_cds;Parent="+query+"-R")
			invert1.append(target+"\tprotein2genome\tmatch_part\t"+str(ee+1)+"\t"+str(es)+"\t.\t"+orient+"\t.\tID="+query+"_matchpart"+";Name="+query+"_matchpart"+";Parent="+query+"-R")
			
			for a in  reverse(invert1):
				print a
print """##FASTA"""