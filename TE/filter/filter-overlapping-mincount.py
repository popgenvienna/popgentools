import sys
import collections
from Modules.SyncTE import SyncTEReader
from optparse import OptionParser, OptionGroup
from optparse import OptionParser,OptionGroup
import copy
import math

#Author: Robert Kofler
parser = OptionParser()
parser.add_option("--te-sync", dest="tesync", help="file containing the sync tes")
parser.add_option("--min-coverage", dest="mincov", help="minimum coverage")
parser.add_option("--ignore-overlap",dest="ign_over",action="store_true",help="Ignore overlapping TE insertions")
parser.add_option("--test",action="store_true", dest="test",help="run the doctest")
(options, args) = parser.parse_args()

mincov=int(options.mincov)
ignoreoverlap=bool(options.ign_over)
tes=SyncTEReader.readfiltered(options.tesync,ignoreoverlap,mincov)

for te in tes:
	print str(te)