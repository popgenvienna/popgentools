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
genomedir=$genomepelementgsnap
outfolder=$2
prefix=$3

echo "using bam ${sam}"
echo "using outfolder ${outfolder}"
echo "using prefix ${prefix}"

# set paths
#pgt use environment variable
#popte use environment variable

s2f="${pgt}/TE/melsim/sam2fastq.py"

# make folders
mkdir -p $outfolder
mkdir -p $outfolder/raw

# let's do it

samtools index $sam
samtools view $sam PPI251 > $outfolder/raw/tmp.sam
python $s2f $outfolder/raw/tmp.sam > $outfolder/raw/tmp.fastq
co1="gsnap -d pele -D ${genomedir} -A sam --novelsplicing=1 -t 4 --quality-protocol illumina ${outfolder}/raw/tmp.fastq | samtools view -Sb - |samtools sort - ${outfolder}/${prefix}.sort"
eval $co1

samtools index $outfolder/${prefix}.sort.bam
rm -rf $outfolder/raw


