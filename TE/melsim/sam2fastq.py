import sys
import random
from optparse import OptionParser, OptionGroup
import collections
import math
import re

parser = OptionParser()
parser.add_option("--sam",dest="sam",help="the input file as sam")
parser.add_option("--count",dest="count",default=False, action="store_true",help="should the thing be counted instead")
(options, args) = parser.parse_args()
count=bool(options.count)

def rc(seq):
        seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
        return "".join([seq_dict[base] for base in reversed(seq)])

c=1
for l in open(options.sam):
        """
        HWUSI-EAS300R:5:6:1167:1537#0	73	PPI251	1	140	74M	*	0	0	CATGATGAAATAACATAAGGTGGTCCCGTCGATAGCCGAAGCTTACCGAAGTATACACTTAAATTCAGTGCACG	a`bbabbab^b_babbbabbab`aaaabaab^babbba\ab`aaaaaa]a`[aaa`]`aaZaa`b^``Ya`XZ\	AS:i:70
        HWUSI-EAS300R:5:9:680:1618#0	99	PPI251	1	11	15S59M	=	167	166	TACGTTTGGTGCGCCCATGATGAAATAACATAAGGTGGTCCCGTCGTAAGCCGAAGCTTACCGAAGTATACACT	a`ab_`aaa^aaaaaaaaaaaa``a`aaaaa`a__Ya_X`_`a_aa^^]]`^aY[__Z[``Z`[[`[\^_^[\a	AS:i:55
        HWUSI-EAS300R:5:11:44:938#0	99	PPI251	1	21	22S52M	=	118	117	TTTGAAGGCGTTTGCTCCACTGCATGATGAAATAACATAAGGTGGTCCCGTCGAAAGCCGAAGCTTACCGAAGT	aaa_aaZ``Y]_a`][][`^_S]]`Y_`_X]]`__^aa__V\X]ZX_]__[Y`X]\]][TXV_]Y]]VW_BBBB	AS:i:52
        HWUSI-EAS300R:5:11:1436:1841#0	163	PPI251	1	13	6S67M1S	=	163	162	AAAGACCATGATGAAATAACATAAGGTGGTCCCGTCGTAAGCCGAAGCTTACCGAAGTATACACTTAAATTCAN	aa\^``_Z_a_]a^\^Z`_W\^XZ`WMVVMT_X]IW^Z^Z_]XUWZZZOXZW]^LLRILXSZZXSWUVPTBBBB	AS:i:63
        HWUSI-EAS300R:5:14:431:796#0	163	PPI251	1	9	22S52M	=	129	128	GGTAAACACAAGTCGACACGACCATGATGAAATAACATAAGGTGGTCCCGTCGTAAGCCGAAGCTTACCGAAGT	ab`aaaaaaaaa_aaaaaaaa`[aabaabaaaa_a_a_a`aaU`aZ]]Y_Y`]U_]_WV^YQ[XX[[\[`^\BB	AS:i:48
        HWUSI-EAS300R:5:14:650:1166#0	163	PPI251	1	28	3S71M	=	91	90	GCGCATGATGAAATAACATAAGGTGGTCCCGTCGATAGCCGAAGCTTACCGAAGTATACACTTAAATTCAGTGC	aabaaabaabaabbabaaa`bab\aa]a`aa^aa_`aa`^a[aa___V[_]V[]Z`^a]^[]\Y\]_\]Z]LY\	AS:i:67
        HWUSI-EAS300R:5:14:1218:856#0	99	PPI251	1	13	8S66M	=	110	109	GGTTCGCGCATGATGAAATAACATAAGGTGGTCCCGTCGGAAGCCGAAGCTTACCGAAGTATACACTTAAATTC	aa_aaa`aa`ab_aa^a]aaa`aaaaaa]aa]aa_a^_a`__`]_a^Y`Z^____aY\]SXW_^_^_]NX[^YB	AS:i:62
        """
        if l.startswith("@"):
                continue
        l=l.rstrip("\n")
        a=l.split("\t")
        readid=a[0]+str(c)
        if count:
                readid=str(c)
                c+=1
        flag=int(a[1])
        seq=a[9]
        qual=a[10]
        if flag & 0x10:
                seq=rc(seq)
                qual="".join([q for q in reversed(qual)])
        print "@{0}".format(readid)
        print seq
        print "+{0}".format(readid)
        print qual
	c+=1
                
        
        
        
