import sys
from optparse import OptionParser,OptionGroup

parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: python igv2manhattan.py -i output.cmh > output
2)	This script transforms your cmh output in a proper input for the R gap library to produce manhattan plots. be aware that this library only accepts numbers as chromosomal names. Thus igv2manhattan transcodes the mel chromosomal names as follows: 
{"2L":"2","2R":"3","3L":"4","3R":"5","4":"6","X":"1"}
	""") 


parser.add_option("-i", "--inp", dest="inp", help="*.sync or cmh output file")


parser.add_option_group(group)
(options, args) = parser.parse_args()

code={"2L":"2","2R":"3","3L":"4","3R":"5","4":"6","X":"1"}

print "CHR\tBP\tP\tgene"
for l in open(options.inp,"r"):
	if "Chromosome" not in l:
		a=l.split()
		if a[0] in code:
			print code[a[0]]+"\t"+a[1]+"\t"+a[-1]+"\tgene"+a[0]+a[1]
