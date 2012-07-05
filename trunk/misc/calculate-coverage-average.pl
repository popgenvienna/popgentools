#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX;


use FindBin qw($RealBin $RealScript);
use lib "$RealBin/../Modules";
use lib "$RealBin/Modules";

#PopGenTools modules
use Pileup;

my $HELP = 0;
my $QUAL_ENCODING = "illumina";
my $MIN_COUNT = 0;
my $MIN_COV = 1;
my $MAX_COV = 1000;
my $MIN_QUAL = 20;

my $CUT_OFF_COVERAGE = 30;
my $MIN_WELL_COVERED_FRACTION = 0;

my $PILEUP_FILE = "";	
my $OUT_FILE = "";
my $CONTIG_LENGTH_FILE="";


GetOptions(
	"pileup=s"=>\$PILEUP_FILE,
	"contigs=s"=>\$CONTIG_LENGTH_FILE,
	"out=s"=>\$OUT_FILE,
	"qual-encoding=s"=>\$QUAL_ENCODING,
	"min-count=i"=>\$MIN_COUNT,
	"min-cov=i"=>\$MIN_COV,
	"max-cov=i"=>\$MAX_COV,
	"min-qual=i"=>\$MIN_QUAL,
	"help"=>\$HELP,
)or die "Invalid arguments, use 'perl calculate-coverage-average.pl --help'.";


pod2usage(-verbose=>3) if $HELP;

pod2usage(-msg=>"Could not find input pileup file.",-verbose=>1) unless defined($PILEUP_FILE);
pod2usage(-msg=>"Could not find contig input file.",-verbose=>1) unless defined($CONTIG_LENGTH_FILE);
pod2usage(-msg=>"Could not find output file.",-verbose=>1) unless defined($OUT_FILE);


print_input_params_append(
$PILEUP_FILE, $OUT_FILE, $CONTIG_LENGTH_FILE, 
$QUAL_ENCODING, $MIN_COV, $MAX_COV, $MIN_QUAL, 
$CUT_OFF_COVERAGE, $MIN_WELL_COVERED_FRACTION);
	
my $ptrContigHash={};

#my $testContigs = "contig1\t12\n".
#				  "contig2\t5\n".
#				  "contig3\t4\n";

#my $testPileup = "contig1\t1\tN\t10\tCCCCCCCCCC\taaaaaaaaaa\n".
#				 "contig1\t2\tN\t15\tCCCCCCCCCCCCCCC\taaaaaaaaaaaaaaa\n".
#				 "contig1\t3\tN\t31\tCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAA\taaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n".
#				 "contig1\t5\tN\t10\tCCCCCCCCCC\taaaaaaaaaa\n".
#				 "contig1\t7\tN\t10\tCCCCCCCCCC\taaaaaaaaaa\n".
#				 "contig1\t10\tN\t10\tCCCCCCCCCC\taaaaaaaaaa\n".
#				 "contig1\t11\tN\t10\tCCCCCCCCCC\taaaaaaaaaa\n".
#				 "contig1\t12\tN\t10\tCCCCCCCCCC\taaaaaaaaaa\n".
#				 "contig2\t1\tN\t10\tCCCCCCCCCC\taaaaaaaaaa\n".
#				 "contig2\t3\tN\t20\tAAAAAAAAAACCCCCCCCCC\taaaaaaaaaaaaaaaaaaaa\n".
#				 "contig2\t5\tN\t10\tCCCCCCCCCC\taaaaaaaaaa\n".
#				 "contig3\t3\tN\t31\tCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\taaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n";

#$CONTIG_LENGTH_FILE = \$testContigs;
#$PILEUP_FILE = \$testPileup;

# load contig names and lengths
open inFileHandle,"<", $CONTIG_LENGTH_FILE;
while (my $line = <inFileHandle>){
	chomp($line);
	
	next if $line =~ m/^#/;
	my ($contigName, $contigLength) = split "\t", $line;
	$ptrContigHash->{$contigName}{length}=$contigLength;
	
}
close inFileHandle;

##parsing pileup file
my $pileupParser=get_pileup_parser($QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL);
open inFileHandle, "<", $PILEUP_FILE;
while (my $line = <inFileHandle>){
	chomp($line);
	my $parsedLine = $pileupParser->($line);
	my $contigName = $parsedLine->{chr};
	
	if ($parsedLine->{eucov}<= $CUT_OFF_COVERAGE){
		$ptrContigHash->{$contigName}{sumCov}+=$parsedLine->{eucov};
	}else{
		$ptrContigHash->{$contigName}{posCoverageOverCutOff}+=1;
		$ptrContigHash->{$contigName}{sumCovOverCovered}+=$parsedLine->{eucov};
	}
	
}
close inFileHandle;

# print out avg coverage for each contig 
open outFileHandle, ">", $OUT_FILE;
print outFileHandle "#contig_name\tcontig_length\taverage_coverage_without_overcovered_positions\taverage_coverage_over_all_positions\tover-covered_fraction\n";

#foreach my contig:
foreach my $contigName (keys %$ptrContigHash){
	my $contigLength = $ptrContigHash->{$contigName}{length};
	my $posCoverageOverCutOff = $ptrContigHash->{$contigName}{posCoverageOverCutOff};

	my $sumCov = $ptrContigHash->{$contigName}{sumCov};
	my $sumCovOverCovered = $ptrContigHash->{$contigName}{sumCovOverCovered};
	
	if (!defined($posCoverageOverCutOff)){
		$posCoverageOverCutOff = 0;
	}

	if (!defined($sumCov)){
		$sumCov=0;
	}

	if (!defined($sumCovOverCovered)){
		$sumCovOverCovered=0;
	}
	
	
	if ( ($contigLength-$posCoverageOverCutOff)/$contigLength >= $MIN_WELL_COVERED_FRACTION){
		my $averageCoverage;
		if ($contigLength == $posCoverageOverCutOff){
			$averageCoverage=0;
		}else{
			$averageCoverage = $sumCov/($contigLength-$posCoverageOverCutOff);
		}
		my $averageCoverageOverAll = ($sumCov+$sumCovOverCovered)/$contigLength;
		print outFileHandle $contigName."\t".$contigLength."\t".$averageCoverage."\t".$averageCoverageOverAll."\t".$posCoverageOverCutOff/$contigLength."\n";  	
	}
}
print "script finished at: ",POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime);
close outFileHandle;



# use right after reading parameters -- you want to have time stamp from the begining of run not from the end ;)
sub print_input_params_append{
	my ($PILEUP_FILE, $OUT_FILE, $CONTIG_LENGTH_FILE, $QUAL_ENCODING, $MIN_COV, $MAX_COV, $MIN_QUAL, $CUT_OFF_COVERAGE, $MIN_WELL_COVERED_FRACTION)=@_;
	
	
	my $out_file_name = "calculate-coverage-average.params";
	my ($inName, $path)= fileparse($PILEUP_FILE);
	$out_file_name = $path.$out_file_name;
	
	open fileHandle, ">>", $out_file_name; 
	print "Info about runned parameters is stored in file $out_file_name and next info will be appended to the file if the script will run again.\n";
	print fileHandle "script $RealBin/$RealScript runned at " , POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime) , " with the following parameter setting:\n";
	print fileHandle "input pileup file: $PILEUP_FILE\n";
	print fileHandle "input contig file (name and length): $CONTIG_LENGTH_FILE";
	print fileHandle "output (contig name, average coverage, overcovered fraction of contig): $OUT_FILE\n";
	print fileHandle "qual-encoding: $QUAL_ENCODING,\n";
	print fileHandle "min-cov: $MIN_COV\n";
	print fileHandle "max-cov: $MAX_COV\n";
	print fileHandle "min-qual: $MIN_QUAL\n";
	print fileHandle "cut-off value for over-coverage: $CUT_OFF_COVERAGE\n";
	print fileHandle "minimal contig fraction that shouldn't be over-coveraged, otherwise the contig will not be in output file: $MIN_WELL_COVERED_FRACTION\n";
	print fileHandle "\n";
	close fileHandle;
}

=head1 NAME

calculate-coverage-average.pl

=head1 SYNOPSIS

perl calculate-coverage-average.pl --pileup FILE_NAME1 --contigs FILE_NAME2 --out FILE_NAME3 

--qual-encoding illumina --min-count 0 --min-cov 1 --max-cov 2000 --min-qual 20 

--coverage-cut-off 30 --min-cov-fraction 0

=head1 DETAILS

The script calculates average coverage for each of contigs from input contigs file. For troublefree run of the script you need to provide names 
of two input files and an output. If you do not set any other parameters the script will run with default settings. 

=head2 OPTIONS

The required options:

=over 10

=item B<--pileup> - An input pileup file.

=item B<--contigs> - An input tab separated file with two items in each row. The first item in a row is a contig name and the second is real length of contig. 
					 You cen obtain this data for example from a header of sam file that corresponds to the pileup file.
					  	
=item B<--out> - An output tab separated file that has three items in each row. The first item is a contig name, the second is an average coverage of the contig 
				and the third item is a fraction of the contig that is over-covered (if it is less than 1 - min-cov-fraction, otherwise the contig will not be 
				in output file). 
				
=back

The other options (with default values):

=over 10

=item B<--qual-encoding> default=illumina.  

=item B<--min-count> default=0. Minimum count of nucleotides of at least two types needed for considering the position as a SNP.

=item B<--min-cov> default=1. Minimum coverage of the any position to be considered as a SNP. 

=item B<--max-cov> default=2000. Maximum coverage of the any position to be considered as a SNP. This should be set at least to the real max coverage in the pileup for this script.

=item B<--min-qual> default=20.


=item B<--coverage-cut-off> default=30. 

=item B<--min-coverage-fraction> default=0.

=back

