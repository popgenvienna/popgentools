#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

use FindBin qw($RealBin);
use lib "$RealBin/../Modules";

use HkaTestModule;

my $GTF_FILE="../../data/dmel-r5-2L.gtf";
my $PILEUP_FILE="../../data/dmel-2L-q20.pileup";
my $MAUVE_PARSED_FILE="../../data/dmel-r5.33+dsim-r1.3-2L.joined";

my $OUT_FILE;#="/Users/tinavisnovska/data/HKAareas.first";

my $QUAL_ENCODING="illumina";
my $MIN_COUNT=2;
my $MIN_COV=10;
my $MAX_COV=10000;
my $MIN_QUAL=20;
my $POOL_SIZE=500;

my $MIN_LENGTH=300;	

my $GENES_COUNT = 50;

##
## loading precalculated measures and overlapping genes for further stat. analysis
##
	
my $VERSION = "v7";	
my $inMeasures = "../../data/$VERSION/Dmel-r5-2L-$VERSION-measures-for-genes";
my $inOverlaps = "../../data/$VERSION/Dmel-r5-2L-$VERSION-overlapping-genes";	
my $outLog="../../data/$VERSION/Dmel-r5-2L-$VERSION.log";


my $ptrMeasures;
my $ptrOverlappingGeneIDs;

($ptrMeasures, $ptrOverlappingGeneIDs) = load_previously_calculated_measures_and_overlaps($inMeasures, $inOverlaps);

my $ptrRemovedGenes;

($ptrMeasures, $ptrRemovedGenes) = handling_with_overlapping_and_short_genes($ptrMeasures, $ptrOverlappingGeneIDs, $MIN_LENGTH);

my $ptrRandomlySelectedGenes=[];

$ptrRandomlySelectedGenes = randomly_choose_genes($ptrMeasures, $GENES_COUNT);

my $T=0;
$T = calculate_T($ptrMeasures,$ptrRandomlySelectedGenes);

calculate_D_expectations_and_variance($ptrMeasures,$ptrRandomlySelectedGenes,$T);

my $chi=0;

$chi=calculate_chi_square_statistic($ptrMeasures, $ptrRandomlySelectedGenes);

print $chi;

my $out = "../../data/$VERSION/dataForStatistic-$VERSION-1";

open my $filehandle, ">", $out;
print_data_used_in_chi_square_statistics($ptrMeasures, $ptrRandomlySelectedGenes, $T, $filehandle);
close $filehandle;


#randomly select $GENES_COUNT genes
#calculate T
#calculate chi square statistics

### open log for appending, print content of $ptrRemovedGenes



##
## calculating measures from scretch (GTF, PILEUP, MAUVE) files needed
##	

##              BEGIN
#my ($ptrMeasures, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs) = load_GPM_collect_genes_calculate_measures_and_overlaps(
#	$GTF_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
#	$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
#	$POOL_SIZE
#);
	
#my $VERSION = "v7";	
#my $outMeasures = "../../data/Dmel-r5-2L-$VERSION-measures-for-genes";
#my $outOverlaps = "../../data/Dmel-r5-2L-$VERSION-overlapping-genes";	
#my $outLog="../../data/Dmel-r5-2L-$VERSION.log";

#open my $outM, ">", $outMeasures;
#open my $outO, ">", $outOverlaps;
#open my $outL, ">", $outLog;

#print_all($ptrMeasures, $ptrSplittedGenes,$ptrOverlaps, $ptrOverlappingGeneIDs, $outM, $outO, $outL);

#close $outL;
#close $outO;
#close $outM;

##           END