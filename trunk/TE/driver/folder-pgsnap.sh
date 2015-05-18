#!/bin/zsh
# read parameters
if [ $# -ne 2 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 infolder outfolder"
exit 2
fi

source ~/.zshrc
infolder=$1
outfolder=$2



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
	torun="zsh ${pgt}/TE/driver/driver-p-gsnap.py ${f} ${subdir} ${prefix}"
	echo $torun
	eval $torun
done
