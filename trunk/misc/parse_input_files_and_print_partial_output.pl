#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use lib "$RealBin/../Modules/External";

use HkaTestModule qw(
load_GPM_collect_genes_calculate_measures_and_overlaps
print_first_round_of_data
);

my $GTF_FILE="../../data/dmel-r5-2L.gtf";
my $PILEUP_FILE="../../data/dmel-2L-q20.pileup";
my $MAUVE_PARSED_FILE="../../data/dmel-r5.33+dsim-r1.3-2L.joined";
my $WITHDRAWN_FILE="../../data/Dmel-withdrawn.txt";

my $OUT_FILE;

my $QUAL_ENCODING="illumina";
my $MIN_COUNT=2;
my $MIN_COV=10;
my $MAX_COV=10000;
my $MIN_QUAL=20;
my $POOL_SIZE=500;

my $MIN_LENGTH=300;	


##
## calculating measures from scretch (GTF, PILEUP, MAUVE) files needed
##	

##              BEGIN
my ($ptrGeneIDs, $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf) = load_GPM_collect_genes_calculate_measures_and_overlaps(
	$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
	$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
	$POOL_SIZE
);
	
my $VERSION = "v12";	
my $outMeasures = "../../data/Dmel-r5-2L-$VERSION-measures-for-genes";
my $outOverlaps = "../../data/Dmel-r5-2L-$VERSION-overlapping-genes";	
my $outLog="../../data/Dmel-r5-2L-$VERSION.log";

open my $outM, ">", $outMeasures;
open my $outO, ">", $outOverlaps;
open my $outL, ">", $outLog;

print_first_round_of_data($ptrGeneIDs, $ptrGeneData, $ptrSplittedGenes,$ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf, $outM, $outO, $outL);

close $outL;
close $outO;
close $outM;

##           END