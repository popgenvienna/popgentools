#!/usr/bin/env python


import sys
from optparse import OptionParser, OptionGroup
import collections
import re
parser = OptionParser()
parser.add_option("--input", dest="input", help="the input file in embl")
parser.add_option("--output-fasta",dest="outfasta",help="output file for the TE sequences in fasta")
parser.add_option("--output-hierarchy",dest="outhier",help="output file for the TE hierarchy file")
(options, args) = parser.parse_args()

def extract_FlyBaseID(toextract):
        # DR   FLYBASE; FBgn0000005; 297.
        tmp=re.split(r"\s+",toextract)
        fbid=tmp[2].rstrip(";")
        teid=tmp[3].rstrip(".")
        # famid  Dsil\Loa
        famid=re.sub(r"D[^\\]+\\","",teid)
        #teid=re.sub(r"\\",";",teid)
        return (fbid,teid,famid)
        
def extract_ID(toextract):
        # ID   DM23420    standard; DNA; INV; 6126 BP.
        tmp=re.split(r"\s+",toextract)
        return tmp[1]
        
def strip_sequence(tostrip):
        #      TACTCGTGCG CATAATAATG CCTTAAAATT TGTTGATGAC ATTCACTCAG TGCAAACAAT       600
        tostrip=re.sub(r"\d+$","",tostrip)
        tostrip=re.sub(r"\s+","",tostrip)
        return tostrip
        
def read_hierarchy(rawhier):
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

def print_hier(outfile,lowHier,highHier):
        # lowhier: ID=>(fbid,teid) #teid=subfamily 
        # highier: teid=>(suborder,order,class)
        ofh=open(outfile,"w")
        ofh.write("insert\tfb\tsubfamily\tfamily\tsuborder\torder\tclass\n")
        testvalue=7
        for k,v in lowHier.items():
                topr=[]
                topr.append(k) # insert
                topr.append(v[0]) # fb
                topr.append(v[1])
                topr.append(v[2])
                key=v[1]
                if(key in highHier):
                        tmp=highHier[key]
                        topr.append(tmp[0])
                        topr.append(tmp[1])
                        topr.append(tmp[2])
                else:
                        raise ValueError("No key found for "+key)
                if len(topr) != testvalue:
                        raise ValueError("length of output is not acceptable, should be "+str(testvalue) + "Output " + "\t".join(topr))
                toprstr="\t".join(topr)
                ofh.write(toprstr+"\n")
        ofh.close()
                
        
        

def print_fasta(outputfile,rawfasta):
        ofh=open(outputfile,"w")
        teh={}
        for fasta in rawfasta:
                fbid=""
                teid=""
                sequence=""
                id=""
                for line in fasta:
                        if(line.startswith("ID")):
                                id=extract_ID(line)
                        elif(line.startswith("DR")):
                                fbid,teid,famid=extract_FlyBaseID(line)
                        elif(line.startswith("     ")):
                                sequence+=strip_sequence(line)+"\n"
                ofh.write(">"+id+"\n")
                ofh.write(sequence)
                if id in teh:
                        raise ValueErro(id +"is present multiple times")
                teh[id]=(fbid,teid,famid)
        return teh
                                



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

lowHier=print_fasta(options.outfasta,rawfasta)
highHier=read_hierarchy(rawhier)
print_hier(options.outhier,lowHier,highHier)
# teh DMTRAM=>(Dmir\TRAM,FBgn0000234)

        
