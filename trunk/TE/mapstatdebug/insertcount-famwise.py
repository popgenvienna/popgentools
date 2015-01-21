from gtfIO import GTFReader,GTFEntry;
import sys
import random
from optparse import OptionParser, OptionGroup
import collections

parser = OptionParser()
parser.add_option("--tes",dest="tes",help="The TE seqs")
(options, args) = parser.parse_args()

teorder=["1360","17.6","1731","297","3S18","412","accord","accord2","aurora-element","baggins","Bari1","Bari2","blood","BS","BS3","BS4","Burdock","Circe","copia","Cr1a","diver","diver2","Dm88","Doc","Doc2-element","Doc3-element","Doc4-element",
	 "F-element","FB","flea","frogger","Fw2","Fw3","G-element","G2","G3","G4","G5","G5A","G6","G7","GATE","gtwin","gypsy","gypsy10","gypsy11","gypsy12","gypsy2","gypsy3","gypsy4","gypsy5",
	 "gypsy6","gypsy7","gypsy8","gypsy9","HB","Helena","HeT-A","HMS-Beagle","HMS-Beagle2","hobo","hopper","hopper2","I-element","Idefix","INE-1","invader1","invader2","invader3","invader4",
	 "invader5","invader6","Ivk","jockey","jockey2","Juan","looper1","Mariner","mariner2","Max-element","McClintock","mdg1","mdg3","micropia","NOF","opus","Osvaldo","P-element","pogo",
	 "Porto1","Q-element","Quasimodo","R1-2","R1A1-element","R2-element","roo","rooA","rover","Rt1a","Rt1b","Rt1c","S-element","S2","springer","Stalker","Stalker2","Stalker3","Stalker4",
	 "Tabor","TAHRE","Tc1","Tc1-2","Tc3","Tirant","Tom1","transib1","transib2","transib3","transib4","Transpac","X-element","ZAM"]



fc=collections.defaultdict(lambda:0)
for l in open(options.tes):
	l=l.rstrip("\n")
	a=l.split("\t")
	# Yhet	43416	F	412	1	LTR	-	-	43296	43316	1	2
	fam=a[3]
	fc[fam]+=1

for fam in teorder:
	count=fc[fam]
	print "{0}\t{1}".format(fam,count)
