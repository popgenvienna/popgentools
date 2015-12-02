#!/bin/zsh
# read parameters
if [ $# -ne 3 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 input.bam outputdir mincount"
exit 2
fi
source ~/.zshrc

bam=$1
of=$2
mincount=$3

#paths
popte2="/Users/robertkofler/dev/popoolation-te2/out/artifacts/popte2_jar/popte2.jar"
hier="/Volumes/Volume_4/analysis/PopTE2/resources/tehier-ml100noS4.fasta"
refg="/Volumes/Volume_4/analysis/PopTE2/refg/chasis1M/chasis1Mtes.fasta"

mkdir $of
mkdir $of/raw


java -jar $popte2 ppileup --bam $bam --map-qual 15 --hier $hier --output $of/raw/pp.gz
java -jar $popte2 identifySignatures --ppileup $of/raw/pp.gz --mode separate --min-count $mincount --output $of/raw/te.signatures
java -jar $popte2 updatestrand --bam $bam --signature $of/raw/te.signatures --output $of/raw/testrand.signatures --hier $hier --map-qual 15 --max-disagreement 0.4
java -jar $popte2 frequency --ppileup $of/raw/pp.gz --signature $of/raw/testrand.signatures --output $of/raw/te.freqsignatures
java -jar $popte2 pairupsignatures --signature $of/raw/te.freqsignatures --ref-genome $refg --hier-file $hier --min-distance -50 --max-distance 300 --output $of/tes.txt

java -jar ~pt2/popte2.jar filterSignatures --input $of/raw/te.freqsignatures --output $of/raw/tefilterd.freqsignatures --max-otherte-count $mincount --max-structvar-count $mincount
java -jar $popte2 pairupsignatures --signature $of/raw/tefilterd.freqsignatures --ref-genome $refg --hier-file $hier --min-distance -50 --max-distance 300 --output $of/tes.filtered.txt

