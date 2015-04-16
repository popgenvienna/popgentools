#!/usr/bin/env python


import sys
from optparse import OptionParser, OptionGroup
import collections
import re
parser = OptionParser()
parser.add_option("--input", dest="input", help="the input file in embl")
parser.add_option("--output-fasta",dest="outfasta",help="output file for the TE sequences in fasta")
parser.add_option("--output-hierarchy",dest="outhier",help="output file for the TE hierarchy file")
parser.add_option("--only-dmeldsim",dest="onlydmeldsim",action="store_true",default=False) # action="store_false", dest="verbose", default=True,
(options, args) = parser.parse_args()


class FastaWriter:
	"""
	Write the content to a fasta file
	"""
	def __init__(self,file,seqleng):
		self.__filename=file
		self.__filehandle=open(file,"w")
		self.__seqleng=seqleng
		
	def write(self,n,s):
		"""
		name, seqeunce
		"""
		sl=self.__seqleng
		fh=self.__filehandle
		fh.write(">"+n+"\n")
		c=0
		while(c<len(s)):
			fh.write(s[c:c+sl]+"\n")
			c+=sl

	def close(self):
		self.__filehandle.close()
	
	@classmethod
	def write_all(cls,file,fastaentries):
		fw=FastaWriter(outputFile)
		for n,s in fastaentries:
			fw.write(n,s)
		fw.close()
		return 1

def read_hierarchy(rawhier):
        """
        FBgn0010103    aurora-element    AB022762          4263bp     ?complete
        """
        resolver={
                "Retroviral":("LTR","LTR","RNA"),
                "non-LTR":("non-LTR","non-LTR","RNA"),
                "SINE-like":("SINE","non-LTR","RNA"),
                "IR-elements:":("IR","TIR","DNA"),
                "MITE":("MITE","TIR","DNA"),
                "Foldback":("Foldback","Foldback","DNA"),
                "Helitron":("Helitron","Helitron","Helitron"),
                "Class":("na","na","na")}
        activeOrder=None
        tehier={}
        for l in rawhier:
                tmp=re.split(r"\s+",l)
                if(tmp[0] in resolver):
                        activeOrder=resolver[tmp[0]]
                elif(tmp[0].startswith("FB")):
                        key=tmp[1]
                        #key=re.sub(r"\\",";",key)
                        tehier[key]=activeOrder
        return tehier

def extract_FlyBaseID(toextract):
        # DR   FLYBASE; FBgn0000005; 297.
        tmp=re.split(r"\s+",toextract)
        fbid=tmp[2].rstrip(";")
        teid=tmp[3].rstrip(".")
        # famid  Dsil\Loa

        return (fbid,teid)
        
def extract_ID(toextract):
        # ID   DM23420    standard; DNA; INV; 6126 BP.
        tmp=re.split(r"\s+",toextract)
        return tmp[1]
        
def strip_sequence(tostrip):
        #      TACTCGTGCG CATAATAATG CCTTAAAATT TGTTGATGAC ATTCACTCAG TGCAAACAAT       600
        tostrip=re.sub(r"\d+$","",tostrip)
        tostrip=re.sub(r"\s+","",tostrip)
        return tostrip


def read_sequences(rawfasta,higherhier):
        """
        ID   DMAURA     standard; DNA; INV; 4263 BP.
        XX
        AC   AB022762;
        XX
        DR   FLYBASE; FBgn0010103; aurora-element.
        XX
        ....
        
        hierarchy
        aurora-element -> [LTR,LTR,RNA]
        """
        entries=[]
        for fasta in rawfasta:
                fbid=""
                teid=""
                sequence=""
                id=""
                for line in fasta:
                        if(line.startswith("ID")):
                                id=extract_ID(line)
                        elif(line.startswith("DR")):
                                fbid,teid=extract_FlyBaseID(line)
                        elif(line.startswith("     ")):
                                sequence+=strip_sequence(line)
                h=[]
                h.append(id)
                h.append(fbid)
                h.append(teid)
                if teid not in higherhier:
                        raise ValueError("key not in higher hierarchyr: "+teid)
                h.extend(higherhier[teid])
                e=Entry(sequence,id,teid,h)
                entries.append(e)
        return entries

def filter_onlyDmelDsim(seqs):
        toret=[]
        for s in seqs:
                if "\\" not in s.family:
                        toret.append(s)
        return toret
                        

class Entry:
        def __init__(self,seq,id,family,hierarchy):
                self.seq=seq
                self.id=id
                self.hierarchy=hierarchy
                self.family=family

def printSequences(seqs,outputfile):
	ofw=FastaWriter(outputfile,60)
	for s in seqs:
		ofw.write(s.id,s.seq)
	ofw.close()            

def printHierarchy(seqs,outputfile):
	ofh=open(outputfile,"w")
	ofh.write("insert\tfb\tfamily\tsuborder\torder\tclass\n")
	for s in seqs:
		tp="\t".join(s.hierarchy)
		ofh.write(tp+"\n")
	ofh.close()

rawhier=[]
rawfasta=[]
flagHierStart=False
flagHierEnd=False
flagSeqStart=False
activeSeq=[]

for line in open(options.input):
        line=line.rstrip("\n")
        if(flagHierEnd==True):
                if(line.startswith("ID")):
                        activeSeq=[line]
                elif(line.startswith("//")):
                        rawfasta.append(activeSeq)
                else:
                        activeSeq.append(line)
                
        elif(flagHierStart==True):
                if(line.startswith("======")):
                        flagHierEnd=True
                else:
                        rawhier.append(line)
        elif(line.startswith("FB gene ID")):
                flagHierStart=True

highHier=read_hierarchy(rawhier) # a hash where the key is a ID like DM23420; the value is a array ["SINE","non-LTR","RNA"]
seqs=read_sequences(rawfasta,highHier)
if(options.onlydmeldsim):
        seqs=filter_onlyDmelDsim(seqs)
        
printHierarchy(seqs,options.outhier)
printSequences(seqs,options.outfasta)
## WRITE ##






        
        


                                

        
