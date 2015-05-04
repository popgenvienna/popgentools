#!/bin/zsh
# read parameters
if [ $# -ne 3 ]
  # "$#" is number of parameters- here we test
  # whether it is not equal to two
   then
   echo "Usage $0 genome.fasta outputdirname prefix"
exit 2
fi
source ~/.zshrc
toa=$1
outfolder=$2
prefix=$3

# set paths
repeatmasker="/Users/robertkofler/programs/RepeatMasker/RepeatMasker"
repeatseq="/Volumes/Volume_4/analysis/male-female-te/resources/te-dmel-dsim-dmau/teseq-dmds-ur.fasta"
repeathier="/Volumes/Volume_4/analysis/male-female-te/resources/te-dmel-dsim-dmau/tehier-dmds-ur.txt"
sciroko="/Users/robertkofler/programs/SciRoKoCo/SciRoKoCo.exe"
pgt="/Users/robertkofler/dev/popgentools"

# make folders
mkdir $outfolder
mkdir $outfolder/rg
mkdir $outfolder/tmp
mkdir $outfolder/results
mkdir $outfolder/popte

# prepare reference genome name
genome="${outfolder}/rg/${prefix}-short.fasta"


# let's do it
cat $toa |awk '{print $1}' > $genome
echo "Start repeatmasking ${genome}"
perl $repeatmasker -gccalc -s -cutoff 200 -no_is -nolow -norna -gff -u -pa 4 -lib $repeatseq $genome
scirout="$outfolder/tmp/${prefix}.td"
echo "Finding microsatellites in ${genome}"
mono $sciroko -mode mmfp -s 12 -p 2 -seedl 8 -seedr 3 -mmao 3 -i $genome -o $scirout
cat $scirout | grep -v '^Seq_Name' | grep -v '^#' | perl -pe 's/^\s+$//'|sed '/^$/d' > ${scirout}.clean
ssr="$outfolder/${prefix}-12-2-8-3.ssr.gff"
cat ${scirout}.clean |awk 'BEGIN{OFS="\t"}{print $1,"SciRoKo","SSR",$4,$5,$7,"+",".","ms"}'> $ssr
echo "Removing microsatellites"
intersectBed -a ${genome}.out.gff -b $ssr -v -f 0.3 > $outfolder/${prefix}.sf.gff
python ${pgt}/TE/annotate/cleanupRepeatMaskerTEannotation.py --mismatch-penalty 0.5 --match-score 1.0 --disable-strand --gtf $outfolder/${prefix}.sf.gff --output $outfolder/${prefix}.resolved.gff

finalannotation="${outfolder}/results/${prefix}.final.gff"
echo "Writing final annotation to ${finalannotation}"
python ${pgt}/TE/annotate/filter-TEannotation-bylength.py --input $outfolder/${prefix}.resolved.gff --min-leng 100 --output $finalannotation

# extract 
python ${pgt}/TE/annotate/mask-repeats.py --input $genome --gtf $finalannotation --output ${outfolder}/results/${prefix}-masked.fasta 
python ${pgt}/TE/annotate/extract-seqhier.py --fasta-ref $genome --gtf-te $finalannotation --hierarchy $repeathier --prefix $prefix --out-fasta ${outfolder}/results/${prefix}-te.fasta  --out-hierarchy ${outfolder}/results/${prefix}-hier.txt

# popte input
cat $repeathier ${outfolder}/results/${prefix}-hier.txt > ${outfolder}/popte/${prefix}-combined-hier.txt
cat ${outfolder}/results/${prefix}-masked.fasta ${outfolder}/results/${prefix}-te.fasta $repeatseq > ${outfolder}/popte/${prefix}-tecombined.fasta






