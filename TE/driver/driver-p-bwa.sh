#!/bin/zsh
# read parameters
if [ $# -ne 3 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 input.bam outputdirname prefix"
exit 2
fi


source ~/.zshrc
set -o shwordsplit

sam=$1
genomefile=$genomepelementbwa
outfolder=$2
prefix=$3

echo "using bam ${sam}"
echo "using outfolder ${outfolder}"
echo "using prefix ${prefix}"

# set paths
#pgt use environment variable
#popte use environment variable

s2f="${pgt}/TE/melsim/sam2fastq.py"
c2ss="${pgt}/TE/melsim/chimera2splicesam.py

# make folders
mkdir -p $outfolder
mkdir -p $outfolder/raw

# let's do it

samtools index $sam
samtools view $sam PPI251 > $outfolder/raw/tmp.sam
python $s2f --count --sam $outfolder/raw/tmp.sam > $outfolder/raw/tmp.fastq
co1="bwa bwasw -t 8 $genomefile $outfolder/raw/tmp.fastq> $outfolder/raw/tmp.sam"
eval co1
python $c2ss -sam $outfolder/raw/tmp.sam |samtools view -T $genomefile -Sb - |samtools sort - ${outfolder}/${prefix}.sort

samtools index $outfolder/${prefix}.sort.bam
rm -rf $outfolder/raw


