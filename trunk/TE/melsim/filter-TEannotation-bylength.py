from gtfIO import GTFReader,GTFWriter;
import sys
import random
from optparse import OptionParser, OptionGroup
import collections

parser = OptionParser()
parser.add_option("--input",dest="input",help="A gtf file containing the RepeatMasked gtf annotation")
parser.add_option("--min-leng",dest="minleng",help="minimum length")
parser.add_option("--output",dest="output",help="A gtf output file")
(options, args) = parser.parse_args()

minleng=int(options.minleng)
w=GTFWriter(options.output)
for e in GTFReader(options.input):
	leng=(e.end-e.start)+1
	if(leng>=minleng):
		w.write(e)
w.close()







