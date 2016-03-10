#!/bin/zsh
# read parameters
if [ $# -ne 4 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 input.bam outputdir mincount expectedstat"
exit 2
fi
source ~/.zshrc

bam=$1
of=$2
mincount=$3
expected=$4

#paths
popte2="/Users/robertkofler/dev/popoolation-te2/out/artifacts/popte2_jar/popte2.jar"
hier="/Volumes/Volume_4/analysis/PopTE2/resources/tehier-ml100noS4.fasta"
refg="/Volumes/Volume_4/analysis/PopTE2/refg/chasis1M/chasis1Mtes.fasta"
scripts="/Volumes/Volume_4/analysis/PopTE2/analysis/2015-11-simulations/scripts"

mkdir $of
mkdir $of/raw


java -jar $popte2 ppileup --bam $bam --map-qual 15 --hier $hier --output $of/raw/pp.gz
java -jar $popte2 identifySignatures --ppileup $of/raw/pp.gz --mode separate --min-count $mincount --output $of/raw/te.signatures
java -jar $popte2 updatestrand --bam $bam --signature $of/raw/te.signatures --output $of/raw/testrand.signatures --hier $hier --map-qual 15 --max-disagreement 0.4
java -jar $popte2 frequency --ppileup $of/raw/pp.gz --signature $of/raw/testrand.signatures --output $of/raw/te.freqsignatures
java -jar $popte2 pairupsignatures --signature $of/raw/te.freqsignatures --ref-genome $refg --hier-file $hier --min-distance -200 --max-distance 300 --output $of/tes.txt

java -jar ~pt2/popte2.jar filterSignatures --input $of/raw/te.freqsignatures --output $of/raw/tefilterd.freqsignatures --max-otherte-count $mincount --max-structvar-count $mincount
java -jar $popte2 pairupsignatures --signature $of/raw/tefilterd.freqsignatures --ref-genome $refg --hier-file $hier --min-distance -200 --max-distance 300 --output $of/tes.filtered.txt


python $scripts/compare-expected-observed.py --expected $expected --observed $of/tes.txt --range 300 > $of/compare.txt
python $scripts/summary-expobs.py $of/compare.txt > $of/summary.txt


python $scripts/compare-expected-observed.py --expected $expected --observed $of/tes.filtered.txt --range 300 > $of/compare.filtered.txt
python $scripts/summary-expobs.py $of/compare.filtered.txt > $of/summary.filtered.txt

