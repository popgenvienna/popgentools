#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
DIR=`pwd`
BAM_NAME=`basename $1 .bam`
[ -e ${BAM_NAME}.bam ] || exit 1
$SAMTOOLS view -h  ${BAM_NAME}.bam | perl -F'\t' -lane 'if (exists($F[10])) {@q = split(//,$F[10]); $F[10] = join("",map {chr(ord($_)-31)} @q)}; print join("\t",@F);'  | $SAMTOOLS  view -bSh - > ${BAM_NAME}_sanger.bam
