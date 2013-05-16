import pysam 
import sys
import argparse
parser = argparse.ArgumentParser(description='converts the quality encoding in bam files from Illumina to Sanger. does no checking, so use with care. needs pysam.') 
parser.add_argument("--in", dest="infile", help="bam file to read", required=True)
parser.add_argument("--out", dest="outfile", help="bam file to output (default: (INFILE wo .bam)_sanger.bam)",default="")
args = parser.parse_args()
infile = vars(args)['infile']
outfile = vars(args)['outfile']

if (outfile==""):
    # strip .bam from right
    outfile=infile[0:infile.rfind(".bam")]
    outfile+="_sanger.bam"
# translation table to convert illumina to sanger
trtbl=''.join([ chr(i-31) if (i-31>0) else chr(0) for i in range(256)])
samfile=pysam.Samfile(infile,"rb")
out=pysam.Samfile(outfile,"wb",template=samfile)
# read lines and convert bam
for l in samfile.fetch(until_eof=True):
    l.qual=l.qual.translate(trtbl)
    out.write(l)
out.close()
samfile.close()
