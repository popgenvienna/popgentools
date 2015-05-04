#!/bin/zsh
# read parameters
if [ $# -ne 3 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 infolder outfolder qualityencoding"
exit 2
fi


source ~/.zshrc
infolder=$1
outfolder=$2
quality=$3


# paths
trimscript="/Volumes/Temp2/Robert/popoolation/basic-pipeline/trim-fastq.pl"

for f in $infolder/*_1.fq.gz
do
 	n=`basename $f`
    	# get rid of file extension
    	n=${n%_1.fq.gz}
	tc="perl $trimscript --min-length 50 --quality-threshold 20  --fastq-type $quality  --input1 $infolder/${n}_1.fq.gz --input2 $infolder/${n}_2.fq.gz --output1 $outfolder/${n}_tr_1.fq.gz --output2 $outfolder/${n}_tr_2.fq.gz --outputse $outfolder/${n}_tr_se.fq.gz"
	echo $tc 

done
