#!/bin/zsh
# read parameters
if [ $# -ne 3 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 infolder outfolder reference"
exit 2
fi


source ~/.zshrc
infolder=$1
outfolder=$2
reference=$3
samro=""


# folders
mkdir $outfolder
mkdir $outfolder/raw

# paths
for f in $infolder/*_1.fq.gz
do
 	n=`basename $f`
    	# get rid of file extension
    	n=${n%_1.fq.gz}
	m1="bwa bwasw -t 12 $reference  $infolder/${n}_1.fq.gz > $outfolder/raw/${n}_1.sam  2> /dev/null"
	m2="bwa bwasw -t 12 $reference  $infolder/${n}_2.fq.gz > $outfolder/raw/${n}_2.sam  2> /dev/null"
	m3="perl /Volumes/Temp2/Robert/popoolationte/samro.pl --fq1 $infolder/${n}_1.fq.gz --fq2 $infolder/${n}_2.fq.gz --sam1 $outfolder/raw/${n}_1.sam --sam2 --output /dev/stdout | samtools view -Sb - | samtools sort -m 4G - $outfolder/${n}.sort"
	echo $m1
	echo $m2
	echo $m3 
done
