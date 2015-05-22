#!/bin/zsh
# read parameters
if [ $# -ne 4 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 input.bam outputdirname prefix quality"
exit 2
fi


set -o shwordsplit


sam=$1
genomedir=$genomepelementgsnap
outfolder=$2
prefix=$3
quality=$4

echo "using bam ${sam}"
echo "using outfolder ${outfolder}"
echo "using prefix ${prefix}"
echo "using quality ${quality}"
echo "using debug ${debug}"

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
python $s2f --sam $outfolder/raw/tmp.sam > $outfolder/raw/tmp.fastq
co1="gsnap -d pele -D ${genomedir} -A sam --novelsplicing=1 -t 4 --quality-protocol ${quality} ${outfolder}/raw/tmp.fastq > $outfolder/raw/tmpgsnap.sam" 
if [ $debug >0 ]
then
echo $co1
fi 
eval $co1


samtools view -Sb  $outfolder/raw/tmpgsnap.sam |samtools sort - ${outfolder}/${prefix}.sort

samtools index $outfolder/${prefix}.sort.bam

if [ $debug -lt 1 ] 
then 
rm -rf $outfolder/raw
fi

