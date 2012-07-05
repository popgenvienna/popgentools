#!/usr/bin/perl

use warnings;
use strict;

use FindBin qw($RealBin $RealScript);
use lib "$RealBin/../Modules";

use Pileup;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX;


my $QUAL_ENCODING="illumina";
my $MIN_COUNT = 0;
my $MIN_COV = 4;
my $MAX_COV = 30;
my $MIN_QUAL = 20;

#testing data:
#my $pileupString="2L\t1\tN\t8\tCcGgcCcG\taaaaaaaa\n".
#				  "2L\t2\tN\t8\tCcGgcCcG\taaaaaaaa\n".
#				  "2L\t3\tN\t8\tCcGgcCcG\taaaaaaaa\n".
#				  "2L\t4\tN\t8\tNNccNNnN\taaaaaaaa\n".
#				  "2L\t5\tN\t8\tAaGgcTTC\taaaaaaaa\n".
#				  "2R\t1\tN\t8\tCcGgcCcG\taaaaaaaa\n".
#				  "2R\t2\tN\t8\tCcGgcCcG\taaaaaaaa\n".
#				  "2R\t3\tN\t8\tCcGgcCcG\taaaaaaaa\n".
#				  "2R\t4\tN\t8\tCcGgcCcG\taaaaaaaa\n".
#				  "2R\t5\tN\t8\tCcGgcCcG\taaaaaaaa\n";

my $PILEUP_FILE="";#\$pileupString
my $ML_RHO="";
my $help;

#read options here:

GetOptions(
	"pileup=s"=>\$PILEUP_FILE,
	"mlRho=s"=>\$ML_RHO,
	"qual-encoding=s"=>\$QUAL_ENCODING,	
	"min-cov=i"=>\$MIN_COV,
	"max-cov=i"=>\$MAX_COV,
	"min-qual=i"=>\$MIN_QUAL,
  	"help"=>\$help
) or die "Invalid arguments, use 'perl pileup2mlRho.pl --help'.";


pod2usage(-verbose=>3) if $help;

pod2usage(-msg=>"Could not find input pileup file.",-verbose=>1) unless defined($PILEUP_FILE);


my $outHandle;
#open output
#order of bases in mlRho input format: ACGT
my $outString;
if ($ML_RHO eq ""){
	open $outHandle, ">-";
	$outString="STDOUT";
}else{
	open $outHandle, ">", $ML_RHO or die "Could not open output file $ML_RHO";
	$outString=$ML_RHO;
};


print_input_params_append($PILEUP_FILE, $outString, $QUAL_ENCODING, $MIN_COV, $MAX_COV, $MIN_QUAL);


my $prevContig="";
my $pileupParser = get_pileup_parser($QUAL_ENCODING, $MIN_COUNT, 0, $MAX_COV, $MIN_QUAL);


open PILEUPfile, "<", $PILEUP_FILE or die "Could not open pileup file $PILEUP_FILE"; 
while (my $line = <PILEUPfile>){
	chomp($line);
	my $parsedLine = $pileupParser->($line);
	my $actualContig = $parsedLine->{chr};
	
	my $position = $parsedLine->{pos};
	my $A = $parsedLine->{A};
	my $C = $parsedLine->{C};
	my $G = $parsedLine->{G};
	my $T = $parsedLine->{T}; 
	
	my $cov = $parsedLine->{eucov};
	
	next unless (($MIN_COV<=$cov)and($cov<=$MAX_COV));

 	if($prevContig eq $actualContig){
	#continue writing output for the same contig	
		print $outHandle $position, "\t", $A, "\t", $C, "\t", $G, "\t", $T, "\n";  
	}else{
	#write a header of a new contig and the first line	
		print $outHandle ">", $actualContig, "\n"; 
		print $outHandle $position, "\t", $A, "\t", $C, "\t", $G, "\t", $T, "\n";
		#update value of prevContig
		$prevContig=$actualContig;
	}
}
close PILEUPfile;
close $outHandle;


# use right after reading parameters -- you want to have time stamp from the begining of run not from the end ;)
sub print_input_params_append{
	my ($IN_FILE, $OUT_FILE, $QUAL_ENCODING, $MIN_COV, $MAX_COV, $MIN_QUAL)=@_;
	
	my $out_file_name = "pileup2mlRho.params";
	my ($inName, $path)= fileparse($IN_FILE);
	$out_file_name = $path.$out_file_name;
	
	open fileHandle, ">>", $out_file_name; 
	print "Info about runned parameters is stored in file $out_file_name and next info will be appended to the file if the script will run again.\n";
	print fileHandle "script $RealBin/$RealScript runned at " , POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime) , " with the following parameter setting:\n";
	print fileHandle "input pileup file: $IN_FILE\n";
	print fileHandle "output (file in an input format for mlRho): $OUT_FILE\n";
	print fileHandle "qual-encoding: $QUAL_ENCODING,\n";
	print fileHandle "min-cov: $MIN_COV\n";
	print fileHandle "max-cov: $MAX_COV\n";
	print fileHandle "min-qual: $MIN_QUAL\n\n";
	close fileHandle;
}

=head1 NAME

pileup2mlRho.pl

=head1 SYNOPSIS

perl pileup2mlRho.pl --pileup _file_name_ --mlRho _file_name_ --qual-encoding illumina --min-cov 2 --max-cov 100 --min-qual 20 

The script transforms pileup data from input to a file format that is required by program mlRho. The only required parameter of this script is pileup file name, all other parameters are optional.
If --mlRho parameter value is not specified, output data will be printed to standard output. B<Default> values of other parameters follow:

qual-encoding: "illumina"

min-cov: 4  -- Minimum coverage per site. It means that sum of counts of nucleotide bases in output will be at least min-cov  

max-cov: 30 -- Maximum coverage per site. See min-cov.

min-qual: 20
	