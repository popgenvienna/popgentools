#!/usr/bin/perl

use warnings;
use strict;

use FindBin qw($RealBin $RealScript);
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX;

use lib "$RealBin/Modules/";
use lib "$RealBin/../Modules/"; 

use TinaPerl qw(read_calculate_print_windows_and_log
				read_data_candidates_list_calculate_print_windows_log);

my $test=undef; 
my $help=0;

my $SYNC_PILEUP_FILE = undef;
my $CANDIDATES_FILE = undef;
my $MAX_COUNT_FOR_3RD_ALLELE = 4;
my $MIN_COUNT_FOR_SELECTED_ALLELES = 2;
my $WINDOW_SIZE = 10000;

my $OUT = undef;
my $C_FILE = undef;
my $IGNORE = 5;
my $LANE = 1;


GetOptions(
	"sync-pileup=s"=>\$SYNC_PILEUP_FILE,
	"candidates=s"=>\$CANDIDATES_FILE,
	"max-count-3=i"=>\$MAX_COUNT_FOR_3RD_ALLELE,
	"min-count=i"=>\$MIN_COUNT_FOR_SELECTED_ALLELES,
	"window-size=i"=>\$WINDOW_SIZE,
	"out=s"=>\$OUT,
	"recomb-file=s"=>\$C_FILE,
	"ignore=i"=>\$IGNORE,
	"lane=i"=>\$LANE,
	"help"=>\$help,
) or die  "Invalid arguments, use 'perl estimate-diff-in-frequencies.pl --help'.";

pod2usage(-verbose=>3) if $help;

pod2usage(-msg=>"Could not find input sync pileup file.",-verbose=>1) unless defined($SYNC_PILEUP_FILE);

pod2usage(-msg=>"Could not find input file with recombination ratios.",-verbose=>1) unless defined($C_FILE);

if(!defined($OUT)){
	
	my ($inName, $inPath)=fileparse($SYNC_PILEUP_FILE);
	$OUT =  $inPath."/diffInFreq";
}

my $out_params = "estimate-diff-in-frequencies.params";
my ($Name, $path)= fileparse($OUT);
$out_params = $path.$out_params;

print "Info about runned parameters is stored in file $out_params and next info will be appended to the file if the script will run again.\n";
open paramsFileHandle, ">>", $out_params;
print paramsFileHandle "script $RealBin/$RealScript started at ", POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime), " with the following parameter setting:\n";
print paramsFileHandle "\tsync-pileup file: ".$SYNC_PILEUP_FILE."\n";
print paramsFileHandle "\tcandidates file: ".$CANDIDATES_FILE."\n";
print paramsFileHandle "\trecomb-file: ".$C_FILE."\n";
print paramsFileHandle "\tout: ".$OUT."\n";
print paramsFileHandle "\twindow-size: ".$WINDOW_SIZE."\n";
print paramsFileHandle "\tmax-count-3: ".$MAX_COUNT_FOR_3RD_ALLELE."\n";
print paramsFileHandle "\tmin-count: ".$MIN_COUNT_FOR_SELECTED_ALLELES."\n";
print paramsFileHandle "\tignore: ".$IGNORE."\n";
print paramsFileHandle "\tlane: ".$LANE."\n";

if (!defined($CANDIDATES_FILE) && !defined($C_FILE)){
	
}elsif(!defined($CANDIDATES_FILE) && defined($C_FILE)){
	read_calculate_print_windows_and_log($SYNC_PILEUP_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $C_FILE, $OUT, $WINDOW_SIZE, $IGNORE, $LANE);
}elsif(defined($CANDIDATES_FILE) && !defined($C_FILE)){
	
}elsif(defined($CANDIDATES_FILE) && defined($C_FILE)){	
	read_data_candidates_list_calculate_print_windows_log($SYNC_PILEUP_FILE, $CANDIDATES_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $C_FILE, $OUT, $WINDOW_SIZE, $IGNORE, $LANE);
};

print paramsFileHandle "script finished at ", POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime), ".\n\n";

# read SNPs
#
#	end -- SNP position
#	early rep 1 -- 1. pop. col
#	late rep 1 -- 3. pop. col
#   choose two alleles
#	print if more than 2 alleles have count > 5
#	skip if unknown or in/del > 0


=head1 NAME 

estimate-diff-in-frequencies.pl

=head1 SYNOPSIS 

perl estimate-diff-in-frequencies.pl --sync-pileup syncPileupFileName --recomb-file fileWithRecombinationRatios --out outFileName

perl estimate-diff-in-frequencies.pl --sync-pileup syncPileupFileName --candidates candidatesListFileName --recomb-file fileWithRecombinationRatios --out outFileName

=head1 OPTIONS

=head2 Default Setting

=over 10

=item min-count = 2

=item max-count-3 = 4

=item window-size = 10000

=back

