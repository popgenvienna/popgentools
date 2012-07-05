from optparse import OptionParser,OptionGroup
import sys


parser = OptionParser()
group=OptionGroup(parser,"""				D E S C R I P T I O N

1)	usage: explainsamflag.py -f 10
2)      this script prints the explanation of the flag in a human readable format
""") 


parser.add_option("-f", "--flag", dest="flag", help="put flag from SAM file")



parser.add_option_group(group)
(options, args) = parser.parse_args()



"""usage %prog decimal-flag [decimal-flag...]

Explain each flag on the command line in plain English
"""


lstFlags = [
    ("read paired", 0x1),
    ("read mapped in proper pair", 0x2),
    ("read unmapped", 0x4),
    ("mate unmapped", 0x8),
    ("read reverse strand", 0x10),
    ("mate reverse strand", 0x20),
    ("first in pair", 0x40),
    ("second in pair", 0x80),
    ("not primary alignment", 0x100),
    ("read fails platform/vendor quality checks", 0x200),
    ("read is PCR or optical duplicate", 0x400)
    ]
    

def explain_sam_flags(iFlags):
    print "Flag",str(iFlags)+":"
    print "________________"    
    for strFlagName, iMask in lstFlags:
        if iFlags & iMask:
            print strFlagName

explain_sam_flags(int(options.flag))
