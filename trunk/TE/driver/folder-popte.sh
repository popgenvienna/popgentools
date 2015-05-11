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
genome=$3

# settings 
pgt="/Volumes/Temp2/Robert/popgentools"


# folders
mkdir -p $outfolder

# paths
for f in $infolder/*.bam
do
 	n=`basename $f`
    	prefix=${n%.bam}
    	prefix=${prefix%.sort}
	subdir="${outfolder}/${prefix}"
	mkdir -p $subdir	
	# Usage driver-popte.py input.[sam/bam] genome.fasta outputdirname prefix
	torun="zsh ${pgt}/TE/driver/driver-popte.py ${f} ${genome} ${subdir} ${prefix}"
	echo $torun
	eval $torun
done
