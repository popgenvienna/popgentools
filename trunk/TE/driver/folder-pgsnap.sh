#!/bin/zsh
# read parameters
if [ $# -lt 2 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 infolder outfolder quality"
exit 2
fi

infolder=$1
outfolder=$2
q=$3


if [ $3 -gt 0 ]
then
	quality="sanger"
else
	quality="illumina"
fi


# settings 
# pgt from environment
mkdir -p $outfolder

# paths
for f in $infolder/*.bam
do
 	n=`basename $f`
    	prefix=${n%.bam}
    	prefix=${prefix%.sort}
	subdir="${outfolder}/${prefix}/p_gsnap"
	mkdir -p $subdir	
	# Usage driver-popte.py input.[sam/bam] genome.fasta outputdirname prefix
	torun="zsh ${pgt}/TE/driver/driver-p-gsnap.py ${f} ${subdir} ${prefix} ${quality}"
	echo $torun
	eval $torun
done
