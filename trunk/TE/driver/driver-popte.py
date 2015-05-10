#!/bin/zsh
# read parameters
if [ $# -ne 4 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 input.[sam/bam] genome.fasta outputdirname prefix"
exit 2
fi


source ~/.zshrc
set -o shwordsplit

sam=$1
genome=$2
outfolder=$3
prefix=$4
temp1=${genome%.fasta}
hier=${temp1}.hier
bootstrap="25 50 75 100 125 150 200"

echo "using sam ${sam}"
echo "using genome ${genome}"
echo "using hierarchy ${hier}"
echo "using outfolder ${outfolder}"
echo "using prefix ${prefix}"
echo "using bootstrap ${bootstrap}"

# set paths
pgt="/Volumes/Temp2/Robert/popgentools"
popte="/Volumes/Temp2/Robert/popoolationte"

# make folders
mkdir -p $outfolder
mkdir -p $outfolder/tmp
mkdir -p $outfolder/results
mkdir -p $outfolder/bootstrap
mkdir -p $outfolder/stat


# sam/bam
if [[ $sam =~ \.bam$ ]] 
then
echo "converting ${sam} to ${outfolder}/tmp/input.sam"
samtools view $sam > $outfolder/tmp/input.sam
else
ln $sam $outfolder/tmp/input.sam
fi
sam=$outfolder/tmp/input.sam


# let's do it
perl $popte/identify-te-insertsites.pl --input $sam  --te-hierarchy-file $hier --te-hierarchy-level family --min-count 1 --min-map-qual 15 --narrow-range 100 --output $outfolder/tmp/tmp.tedirectional > /dev/null
# N-s
perl $popte/genomic-N-2gtf.pl --input $genome > $outfolder/tmp/polyn.gtf
# crosslink
perl $popte/crosslink-te-sites.pl --directional-insertions $outfolder/tmp/tmp.tedirectional --min-dist 85 --max-dist 400 --output $outfolder/tmp/tmp.teinsertions --poly-n $outfolder/tmp/polyn.gtf --te-hierarchy $hier > /dev/null
# estimate poly

perl $popte/estimate-polymorphism.pl --sam-file mapte/Dsim/Dsim.sort.sam --te-insert-file te/Dsim/raw/Dsim.teinsertions --te-hierarchy-file ~dmau/refg/Dsim/Dsim-combinedte-incM.hier --te-hierarchy-level family --min-map-qual 15 --output $outfolder/tmp/tmp.polytes > /dev/null

# filter major chromosome 
cat $outfolder/tmp/tmp.polytes |awk '$1=="X" || $1=="2L" || $1=="2R" || $1=="3L" || $1=="3R" || $1=="4"' > $outfolder/tmp/tmp.euchr.polytes
perl $pgt/TE/filter/filter-teinserts.pl --te-insertions $outfolder/tmp/tmp.euchr.polytes --min-count 10 --output $outfolder/results/${prefix}.mc10.euchr.polytes
perl $pgt/TE/filter/filter-teinserts.pl --discard-overlapping --te-insertions $outfolder/results/${prefix}.mc10.euchr.polytes --output  $outfolder/results/${prefix}.mc10.nooverlap.euchr.polytes 

for b in $bootstrap
do
 python $pgt/TE/melsim/targetcoverage.py --tes $outfolder/results/${prefix}.mc10.euchr.polytes --target-cov $b --min-count 1 > $outfolder/bootstrap/${prefix}.boot${b}.polytes
 perl $pgt/TE/filter/filter-teinserts.pl --discard-overlapping --te-insertions $outfolder/bootstrap/${prefix}.boot${b}.polytes --output $outfolder/bootstrap/${prefix}.boot${b}.nooverlap.polytes
 python $pgt/TE/mapstatdebug/insertcount-famwise.py --tes $outfolder/bootstrap/${prefix}.boot${b}.polytes > $outfolder/stat/${prefix}.boot${b}.inscount.txt
python $pgt/TE/mapstatdebug/insertcount-famwise.py --tes $outfolder/bootstrap/${prefix}.boot${b}.nooverlap.polytes > $outfolder/stat/${prefix}.boot${b}.nooverlap.inscount.txt
done

# statistics
stat1="perl ${pgt}/TE/mapstatdebug/debugstat-mapedreads-familywise.pl --input ${sam} --te-annotation ${hier} > ${outfolder}/stat/${prefix}.mapedreads.txt"
#echo $stat1
eval $stat1
stat2="perl ${pgt}/TE/mapstatdebug/debugstat-presFrag-familywise.pl --input ${sam} --te-annotation ${hier} > ${outfolder}/stat/${prefix}.presfrag.txt"
#echo $stat2
eval $stat2
python $pgt/TE/mapstatdebug/insertcount-famwise.py --tes $outfolder/results/${prefix}.mc10.euchr.polytes > $outfolder/stat/${prefix}.mc10.euchr.inscount.txt
python $pgt/TE/mapstatdebug/insertcount-famwise.py --tes $outfolder/results/${prefix}.mc10.nooverlap.euchr.polytes > $outfolder/stat/${prefix}.mc10.nooverlap.euchr.inscount.txt

#clean 
rm $outfolder/tmp/input.sam


