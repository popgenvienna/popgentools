#!/bin/sh

for d in `ls data/*/*.bin`
do
	DIR=`dirname $d`
	GEN=`basename $DIR`
	
	echo $GEN
	zip -r snpEff_v2_0_2_$GEN.zip $d
done
