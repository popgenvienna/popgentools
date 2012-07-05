#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;
use File::Basename;

use POSIX;
use FindBin qw($RealBin $RealScript);
use lib "$RealBin/../Modules";
use lib "$RealBin/../Modules/External";

use HkaTestModule qw(
load_GPM_collect_genes_calculate_measures_and_overlaps
load_GPM_collect_genes_calculate_measures_and_overlaps_old_theta
load_simulated_GP_collect_genes_calculate_measures_and_overlaps
load_simulated_GP_collect_genes_calculate_measures_and_overlaps_old_theta
print_first_round_of_data


load_previously_calculated_measures_and_overlaps
load_reference_gene_list
handling_with_overlapping_and_short_and_zero_variance_genes
calculate_and_store_chi_square_statistic_for_each_gene
randomly_choose_genes
print_selected_genes_with_all_statistics
print_reference_set_of_genes_with_T_for_the_set
print_log_removed_genes
contains

);


# input raw data
my $GTF_FILE = undef;
my $PILEUP_FILE = undef;
my $MAUVE_PARSED_FILE= undef;
my $WITHDRAWN_FILE= undef;

my $OUT = undef;

# input summary
my $SUMMARY_FILE= undef; 

my $SUMMARY=undef;
my $RAW_DATA=undef;
my $SFRD = undef;

#input both
my $QUAL_ENCODING="illumina";
my $MIN_COUNT=2;
	# min count per base on the site
my $MIN_COV=4;	
	# max count per base on the site
my $MAX_COV=1000;
my $MIN_QUAL=20;
my $POOL_SIZE=500;

my $MIN_LENGTH=700;
my $GENES_COUNT = 100;	
my $REPETITIONS = 50;

my $REFERENCE_GENES_SET_FILE=undef;
	#if defined then take as a reference set	 

my $HELP= 0;
my $TEST= 0;
# testuj GffGtfParser, HKAtestModule

my $VERSION = undef;
my $OLD = 0;
my $NO_MAUVE = 0;

my $GENES_INFO_FILE = undef;
my $OVERLAPS_FILE = undef;

my $SEGREGATING_SITES="avgCorr";

# avgCn
# perSiteCn
# avgCorr
# perSiteCorr

GetOptions(
	"gtf-file=s"=>\$GTF_FILE,
	"pileup-file=s"=>\$PILEUP_FILE,
	"mauve-parsed-file=s"=>\$MAUVE_PARSED_FILE,
	"withdrawn-file=s"=>\$WITHDRAWN_FILE,
	
	"genes-info-file=s"=>\$GENES_INFO_FILE,
	"overlaps-file=s"=>\$OVERLAPS_FILE,
	"reference-genes-set-file=s"=>\$REFERENCE_GENES_SET_FILE,
	"out=s"=>\$OUT,
	
	"min-count=i"=>\$MIN_COUNT,
	
	"repetitions=i"=>\$REPETITIONS,
	"min-length-for-statistic=i"=>\$MIN_LENGTH,
	"reference-set-size=i"=>\$GENES_COUNT,
	
	"version=s"=>\$VERSION,
	#"old"=>\$OLD,
	"no-mauve"=>\$NO_MAUVE,
	
	"segregating-sites=s"=>\$SEGREGATING_SITES,
	
	"statistic-from-summary"=>\$SUMMARY,
	"statistic-from-raw-data"=>\$RAW_DATA,
	"summary-from-raw-data"=>\$SFRD,
	
	"help"=>\$HELP,
);


#TODO: tests

pod2usage(-verbose=>3) if $HELP;

my @input = ($SUMMARY, $RAW_DATA, $SFRD);
#check, that exactly one of the possibilities is choosed:
my $sum=0;
foreach my $element (@input){
	if (defined($element) and ($element == 1)){
		$sum+=1;
	}
}

if (! ($sum == 1) ){
	pod2usage(-msg=>"Specify exactly one of the following options: \n\t--statistic-from-raw-data,\n\t--statistic-from-summary,\n\t--summary-from-raw-data.",-verbose=>1);	
	
}


my $SEGREGATING_SITES_BAD;
my $segregating_sites_changed=0;

if (!($SEGREGATING_SITES eq "perSiteCn") and !($SEGREGATING_SITES eq "avgCn") and !($SEGREGATING_SITES eq "perSiteCorr") and !($SEGREGATING_SITES eq "avgCorr") ){
	$segregating_sites_changed = 1;
	$SEGREGATING_SITES_BAD = $SEGREGATING_SITES;
	$SEGREGATING_SITES = "perSiteCorr";
}


if (defined($SUMMARY) and ($SUMMARY == 1)){
# calculate statistic from summary data 	
	if( ! (defined($GENES_INFO_FILE) and defined($OVERLAPS_FILE) ) ){
		pod2usage(-msg=>"For calculation statistic from previously created summary files, you should specify theese options: \n\t--genes-info-file,\n\t--overlaps-file.",-verbose=>1);
	}else{
		my ($inName, $path)= fileparse($GENES_INFO_FILE);
		my $paramsFileName = $path."HKAtest.params";
		
		if (!defined($OUT)){
			$OUT = $path."HKAout-statFromSummary";
		}
		
		#print input to params
		open my $paramsFileHandle, ">>", $paramsFileName;
		
 		print "Info about runned parameters is stored in file $paramsFileName and next info will be appended to the file if the script will run again.\n";
 		print $paramsFileHandle "Script $RealBin/$RealScript started at " , POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime) , " with the following parameters setting:\n";
 		print $paramsFileHandle "\tStatistic calculated from summary files:\n";
 		print $paramsFileHandle "\t\tgenes-info-file: ", $GENES_INFO_FILE, "\n";
		print $paramsFileHandle "\t\toverlaps-file: ", $OVERLAPS_FILE, "\n";
		print $paramsFileHandle "\tOther setting:\n"; 
	
		my $ptrGeneData={};
		my $ptrOverlappingGeneIDs={};
		my $ptrRemovedGenes={};
		my $ptrGeneIDs=[];
			
		($ptrGeneData, $ptrOverlappingGeneIDs) = load_previously_calculated_measures_and_overlaps($GENES_INFO_FILE, $OVERLAPS_FILE);
		($ptrGeneIDs, $ptrRemovedGenes) = handling_with_overlapping_and_short_and_zero_variance_genes($ptrGeneData, $ptrOverlappingGeneIDs, $MIN_LENGTH);

		calculate_and_print_statistics($ptrGeneData, $ptrGeneIDs, $ptrRemovedGenes, $OUT, $REFERENCE_GENES_SET_FILE, $GENES_COUNT, $REPETITIONS, $paramsFileHandle);
		print $paramsFileHandle "The script finished at " , POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime), "\n\n";
		close $paramsFileHandle;
	}
}elsif(defined($RAW_DATA) and ($RAW_DATA == 1)){
	
	if (!(defined($GTF_FILE) and defined($PILEUP_FILE) and defined($MAUVE_PARSED_FILE) )){	
		pod2usage(-msg=>"For calculation statistics from raw data, you should specify theese options: \n\t--gtf-file,\n\t--pileup-file\n\t--mauve-parsed-file.",-verbose=>1);
	}else{
		my $withdrawnEmpty=0;
		if (!defined($WITHDRAWN_FILE)){my $empty = ""; $WITHDRAWN_FILE=\$empty; $withdrawnEmpty=1;}
		if (!defined($OLD)){$OLD=0}

		my ($inName, $path)= fileparse($GTF_FILE);
		my $paramsFileName = $path."HKAtest.params";
		
		open my $paramsFileHandle, ">>", $paramsFileName;
		
		print "Info about runned parameters is stored in file $paramsFileName and next info will be appended to the file if the script will run again.\n";
 		print $paramsFileHandle "Script $RealBin/$RealScript started at " , POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime) , " with the following parameters setting:\n";
 		print $paramsFileHandle "\tStatistic calculated from raw data files:\n";
		print $paramsFileHandle "\t\tgtf-file: ", $GTF_FILE, "\n";
		print $paramsFileHandle "\t\tpileup-file: ", $PILEUP_FILE, "\n";
		print $paramsFileHandle "\t\tmauve-parsed-file: ", $MAUVE_PARSED_FILE, "\n"; 
		print $paramsFileHandle "\tOther setting:\n";
		if (!$withdrawnEmpty){
			print $paramsFileHandle "\t\twithdrawn-file (list of genes that sohuld be withdrawn): ", $WITHDRAWN_FILE, "\n"; 
		}
		if($segregating_sites_changed and defined($SEGREGATING_SITES_BAD)){
			print $paramsFileHandle "\t\tsegregating-sites: ", $SEGREGATING_SITES, ", changed from: ", $SEGREGATING_SITES_BAD, "\n";
		}else{
			print $paramsFileHandle "\t\tsegregating-sites: ", $SEGREGATING_SITES, "\n";
		}
		print $paramsFileHandle "\t\tqual-encoding: ", $QUAL_ENCODING, "\n";
		print $paramsFileHandle "\t\tmin-count: ", $MIN_COUNT, "\n";
		print $paramsFileHandle "\t\tmin-cov: ", $MIN_COV, "\n";
		print $paramsFileHandle "\t\tmax-cov: ", $MAX_COV, "\n";
		print $paramsFileHandle "\t\tmin-qual: ", $MIN_QUAL, "\n";
		print $paramsFileHandle "\t\tpool-size: ", $POOL_SIZE, "\n";
		print $paramsFileHandle "\t\told (theta): ", $OLD, "\n";
		
#		if ($OLD){
#			if (!defined($OUT)){		
#				$OUT = $path."HKAout-statFromRawData-old-theta";
#			}
#	
#			my ($ptrGeneIDs, 
#				$ptrGeneData, 
#				$ptrSplittedGenes, 
#				$ptrOverlaps, 
#				$ptrOverlappingGeneIDs, 
#				$ptrWithdrawnInGtf
#			) = load_GPM_collect_genes_calculate_measures_and_overlaps_old_theta(
#																		$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
#																		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
#																		$POOL_SIZE, 
#																		$SEGREGATING_SITES 
#																	  );	
#
#			print_summary($ptrGeneIDs,  $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf, $OUT, $paramsFileHandle);
#						
#			my $ptrRemovedGenes={};
#			($ptrGeneIDs, $ptrRemovedGenes) = handling_with_overlapping_and_short_and_zero_variance_genes($ptrGeneData, $ptrOverlappingGeneIDs, $MIN_LENGTH, $ptrGeneIDs);
#			
#			calculate_and_print_statistics($ptrGeneData, $ptrGeneIDs, $ptrRemovedGenes, $OUT, $REFERENCE_GENES_SET_FILE, $GENES_COUNT, $REPETITIONS, $paramsFileHandle);
#		}else{
			if (!defined($OUT)){
				$OUT = $path."HKAout-statFromRawData";
			}
	
			my ($ptrGeneIDs, 
				$ptrGeneData, 
				$ptrSplittedGenes, 
				$ptrOverlaps, 
				$ptrOverlappingGeneIDs, 
				$ptrWithdrawnInGtf
			) = load_GPM_collect_genes_calculate_measures_and_overlaps(
																		$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
																		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
																		$POOL_SIZE, 
																		$SEGREGATING_SITES
																	  );	
			
			print_summary($ptrGeneIDs,  $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf, $OUT, $paramsFileHandle);
			
			my $ptrRemovedGenes={};
			($ptrGeneIDs, $ptrRemovedGenes) = handling_with_overlapping_and_short_and_zero_variance_genes($ptrGeneData, $ptrOverlappingGeneIDs, $MIN_LENGTH, $ptrGeneIDs);
			
			calculate_and_print_statistics($ptrGeneData, $ptrGeneIDs, $ptrRemovedGenes, $OUT, $REFERENCE_GENES_SET_FILE, $GENES_COUNT, $REPETITIONS, $paramsFileHandle);
#		}
		print $paramsFileHandle "The script finished at " , POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime), "\n\n";
		close $paramsFileHandle;
				
	} 
}elsif(defined($SFRD) and ($SFRD==1)){

	if (!(defined($GTF_FILE) and defined($PILEUP_FILE) and defined($MAUVE_PARSED_FILE) )){	
		pod2usage(-msg=>"For calculation summary from raw data, you should specify theese options: \n\t--gtf-file,\n\t--pileup-file\n\t--mauve-parsed-file.",-verbose=>1);
	}else{
		my $withdrawnEmpty=0;
		if (!defined($WITHDRAWN_FILE)){my $empty = ""; $WITHDRAWN_FILE=\$empty; $withdrawnEmpty=1;}
		if (!defined($OLD)){$OLD=0}
	
		my ($inName, $path)= fileparse($GTF_FILE);
		my $paramsFileName = $path."HKAtest.params";
		
		open my $paramsFileHandle, ">>", $paramsFileName;
		
		print "Info about runned parameters is stored in file $paramsFileName and next info will be appended to the file if the script will run again.\n";
 		print $paramsFileHandle "Script $RealBin/$RealScript started at " , POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime) , " with the following parameters setting:\n";
 		print $paramsFileHandle "\tSummary calculated from raw data files:\n";
		print $paramsFileHandle "\t\tgtf-file: ", $GTF_FILE, "\n";
		print $paramsFileHandle "\t\tpileup-file: ", $PILEUP_FILE, "\n";
		print $paramsFileHandle "\t\tmauve-parsed-file: ", $MAUVE_PARSED_FILE, "\n"; 
		print $paramsFileHandle "\tOther setting:\n";
		if (!$withdrawnEmpty){
			print $paramsFileHandle "\t\twithdrawn-file (list of genes that sohuld be withdrawn): ", $WITHDRAWN_FILE; 
		}
		print $paramsFileHandle "\t\tsegregating-sites: ", $SEGREGATING_SITES, "\n";
		print $paramsFileHandle "\t\tqual-encoding: ", $QUAL_ENCODING, "\n";
		print $paramsFileHandle "\t\tmin-count: ", $MIN_COUNT, "\n";
		print $paramsFileHandle "\t\tmin-cov: ", $MIN_COV, "\n";
		print $paramsFileHandle "\t\tmax-cov: ", $MAX_COV, "\n";
		print $paramsFileHandle "\t\tmin-qual: ", $MIN_QUAL, "\n";
		print $paramsFileHandle "\t\tpool-size: ", $POOL_SIZE, "\n";
		print $paramsFileHandle "\t\told (theta): ", $OLD, "\n";
		
#		if ($OLD){
#			if (!defined($OUT)){
#				my ($inName, $path)= fileparse($GTF_FILE);
#				$OUT = $path."HKAout-symmaryFromRawData-old-theta";
#			}
#			
#			my ($ptrGeneIDs, 
#	    	    $ptrGeneData, 
#	    		$ptrSplittedGenes, 
#	    		$ptrOverlaps, 
#	    		$ptrOverlappingGeneIDs, 
#	    		$ptrWithdrawnInGtf
#			) = load_GPM_collect_genes_calculate_measures_and_overlaps_old_theta(
#																		$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
#																		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
#																		$POOL_SIZE,
#																		$SEGREGATING_SITES
#																	  );
#			print_summary($ptrGeneIDs,  $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf, $OUT, $paramsFileHandle);
#		}else{
			if (!defined($OUT)){
				my ($inName, $path)= fileparse($GTF_FILE);
				$OUT = $path."HKAout-symmaryFromRawData";
			}
		
			my ($ptrGeneIDs, 
	    	    $ptrGeneData, 
	    		$ptrSplittedGenes, 
	    		$ptrOverlaps, 
	    		$ptrOverlappingGeneIDs, 
	    		$ptrWithdrawnInGtf
			) = load_GPM_collect_genes_calculate_measures_and_overlaps(
																		$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
																		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
																		$POOL_SIZE,
																		$SEGREGATING_SITES
																	  );
			print_summary($ptrGeneIDs,  $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf, $OUT, $paramsFileHandle);
#		}
		print $paramsFileHandle "The script finished at " , POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime), "\n\n";
		close $paramsFileHandle;
	
	}
}


#load_simulated_GP_collect_genes_calculate_measures_and_overlaps_old_theta($GTF_FILE, $PILEUP_FILE, $QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	$POOL_SIZE);
	
sub print_summary{
	my ($ptrGeneIDs,  $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf, $OUT, $paramsFileHandle)=@_;
	
	my $outMeasures = $OUT.".measuresForGenes";
	my $outOverlaps = $OUT.".overlappingGenes";
	my $outLog = $OUT.".log";
	
	open my $outM, ">", $outMeasures;
	open my $outO, ">", $outOverlaps;
	open my $outL, ">", $outLog;

	print_first_round_of_data($ptrGeneIDs, $ptrGeneData, $ptrSplittedGenes,$ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf, $outM, $outO, $outL);

	close $outL;
	close $outO;
	close $outM;
	
	print $paramsFileHandle "\tOutput files: \n";
	print $paramsFileHandle "\t\tgenes-info-file: ", $outMeasures, "\n";
	print $paramsFileHandle "\t\toverlaps-file: ", $outOverlaps, "\n";
	print $paramsFileHandle "\t\tlog: ", $outLog, "\n"; 
}

sub print_statistic{
	my ($OUT, $ptrGeneIDsFinal, $ptrReferenceSetOfGenes, $ptrGeneData, $ptrRemovedGenes, $paramsFileHandle)=@_;
	
	my $outStat = $OUT.".statistics"; 
	my $outRef = $OUT.".referenceSet";
	my $outLog = $OUT.".log";
	
	#print statistics for all geneIDs
	open my $fileHandle, ">", $outStat;
	print_selected_genes_with_all_statistics($ptrGeneIDsFinal, $ptrGeneData, $fileHandle);
	close $fileHandle;

	#print statistics for ref set of genes
	open $fileHandle, ">", $outRef;
	print_reference_set_of_genes_with_T_for_the_set($ptrReferenceSetOfGenes, $ptrGeneData, $fileHandle);
	close $fileHandle;

	#print log file
	open $fileHandle, ">>", $outLog;
	print_log_removed_genes($ptrRemovedGenes, $fileHandle);
	close $fileHandle;
	
	print $paramsFileHandle "\tOutput files: \n";
	print $paramsFileHandle "\t\tstatistic: ", $outStat, "\n";
	print $paramsFileHandle "\t\treference set: ", $outRef, "\n";
	print $paramsFileHandle "\t\tlog: ", $outLog, "\n";
}

sub calculate_and_print_statistics{
	my ($ptrGeneData, $ptrGeneIDs, $ptrRemovedGenes, $OUT, $REFERENCE_GENES_SET_FILE, $GENES_COUNT, $REPETITIONS, $paramsFileHandle)=@_;
	
	my $ptrReferenceSetOfGenes=[];
	my $ptrGeneIDsFinal=[];
		
	if (defined($REFERENCE_GENES_SET_FILE)){				

		print $paramsFileHandle "\t\treference-genes-set-file: ", $REFERENCE_GENES_SET_FILE, "\n";
		print $paramsFileHandle "\t\tmin-length-for-statistic: ", $MIN_LENGTH, "\n";

		$ptrReferenceSetOfGenes=load_reference_gene_list($REFERENCE_GENES_SET_FILE);
						
		foreach my $geneID (@{$ptrGeneIDs}){
			next if contains($geneID, $ptrReferenceSetOfGenes);
			push @{$ptrGeneIDsFinal}, $geneID;
		}

		calculate_and_store_chi_square_statistic_for_each_gene($ptrGeneData, $ptrReferenceSetOfGenes, $ptrGeneIDsFinal);
		print_statistic($OUT, $ptrGeneIDsFinal, $ptrReferenceSetOfGenes, $ptrGeneData, $ptrRemovedGenes, $paramsFileHandle);			
				 
	}else{
		my $ptrGeneDataBackup=$ptrGeneData;

		print $paramsFileHandle "\t\treference-set-size: ", $GENES_COUNT, "\n";
		print $paramsFileHandle "\t\trepetitions: ", $REPETITIONS, "\n";
		print $paramsFileHandle "\t\tmin-length: ", $MIN_LENGTH, "\n";
			
		for (my $i=1; $i<=$REPETITIONS; $i++){
			$ptrGeneData = $ptrGeneDataBackup;		
			$ptrReferenceSetOfGenes = randomly_choose_genes($ptrGeneData, $ptrGeneIDs, $GENES_COUNT);

			foreach my $geneID (@{$ptrGeneIDs}){
				next if contains($geneID, $ptrReferenceSetOfGenes);
				push @{$ptrGeneIDsFinal}, $geneID;
			}
			
			calculate_and_store_chi_square_statistic_for_each_gene($ptrGeneData, $ptrReferenceSetOfGenes, $ptrGeneIDsFinal);
			my $oout = $OUT.".".$i;	
			print_statistic($oout, $ptrGeneIDsFinal, $ptrReferenceSetOfGenes, $ptrGeneData, $ptrRemovedGenes, $paramsFileHandle);
			$ptrReferenceSetOfGenes=[];
			$ptrGeneIDsFinal=[];				
		}
	}	
}

=head1 NAME

HKAtest.pl

=head1 SYNOPSIS

perl HKAtest.pl --statistic-from-raw-data --gtf-file gtfFileName --pileup-file pileupFileName --mauve-parsed-file mauveParsedFileName --reference-genes-set-file fileName --out fileName

perl HKAtest.pl --statistic-from-raw-data --gtf-file gtfFileName --pileup-file pileupFileName --mauve-parsed-file mauveParsedFileName --withdrawn-file withdrawnList --old --out fileName

perl HKAtest.pl --statistic-from-summary --genes-info-file fileName --overlaps-file fileName --reference-genes-set-file fileName --out fileName

perl HKAtest.pl --statistic-from-summary --genes-info-file fileName --overlaps-file fileName --reference-set-size 100 --min-length-for-statistic 700 --repetitions 50 --out fileName

perl HKAtest.pl --summary-from-raw-data --gtf-file gtfFileName --pileup-file pileupFileName --mauve-parsed-file mauveParsedFileName --reference-genes-set-file fileName --out fileName 




=head1 DEFAULT SETTINGS

=over 10

=item --qual-encoding illumina

=item --min-count 2

=item --min-cov 4 

=item --max-cov 1000

=item --min-qual 20

=item --pool-size 500

=item --min-length-for-statistic 700

=item --reference-set-size 100

=item --repetitions 50

=item --segregating-sites avgCorr

=item 

=item --old 0

=item --withdrawn-file undef

=item --reference-genes-set-file undef

=back

