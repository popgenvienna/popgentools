#!/usr/bin/perl

use warnings;
use strict;

use FindBin qw($RealBin $RealScript);
use lib "$RealBin/../Modules";

#PopGenTools
use TinaPerl;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX;
 
my $help = 0;
my $test="";

# default settings: 

#my $SYNC_VARIANCE_FILE = "/Users/tinavisnovska/data/Pablo/window/BR9BR3mBR4_F15r4mF15r5F23r1_F27r5_F37r4F37r5F37r1_synchronized_4_PI";
#my $OUTPUT_FILE = "/Users/tinavisnovska/data/Pablo/window/out2";
# change!!!!!
#
my $SYNC_VARIANCE_FILE = undef;
my $OUTPUT_FILE = undef;

my $WINDOW_SIZE = 1000;
my $STEP_SIZE = 1000;
my $MIN_COV_FRACTION = 0.6;
my $MIN_LENGTH_FRACTION = 1;
my $ALL_POP_NONZERO = 0;
my $PRINT_ZERO_LINES = 0; 
 
GetOptions(
	"sync-variance=s"=>\$SYNC_VARIANCE_FILE,
	"output=s"=>\$OUTPUT_FILE,
	"window-size=i"=>\$WINDOW_SIZE,
	"step-size=i"=>\$STEP_SIZE,
	"min-cov-fraction=i"=>\$MIN_COV_FRACTION,
	"min-length-fraction=i"=>\$MIN_LENGTH_FRACTION,
	"all-pop-nonzero"=>\$ALL_POP_NONZERO,
	"print-zero-lines"=>\$PRINT_ZERO_LINES,
	"test=s"=>\$test,
  	"help"=>\$help
) or die "Invalid arguments, use 'perl variance-for-sync.pl --help'.";


if ($test eq "prove"){
	my $testOut;
	$testOut = `prove $RealBin/../Modules/Test/variance-for-sync-sliding-window.t`;
	print $testOut, "\n";
}elsif ($test eq "perl"){
	my $testOut;
	$testOut = `perl $RealBin/../Modules/Test/variance-for-sync-sliding-window.t`;	
	print $testOut, "\n";
}else{
	pod2usage(-verbose=>3) if $help;

	pod2usage(-msg=>"Could not find input sync variance file.",-verbose=>1) unless defined($SYNC_VARIANCE_FILE);

	pod2usage(-msg=>"Could not find name of output file, it is necessary parameter, please specify it.", -verbose=>1) unless defined($OUTPUT_FILE);

	print_input_params_append($SYNC_VARIANCE_FILE, $OUTPUT_FILE, $WINDOW_SIZE, $STEP_SIZE, $MIN_COV_FRACTION, $MIN_LENGTH_FRACTION, $ALL_POP_NONZERO);

	average_variance_sliding_window_sync($SYNC_VARIANCE_FILE, $OUTPUT_FILE, $WINDOW_SIZE, $STEP_SIZE, $MIN_COV_FRACTION, $MIN_LENGTH_FRACTION, $PRINT_ZERO_LINES, $ALL_POP_NONZERO);
}


sub print_input_params_append{
	my ($IN_FILE, $OUTPUT_FILE, $WINDOW_SIZE, $STEP_SIZE, $MIN_COV_FRACTION, $MIN_LENGTH_FRACTION, $ALL_POP_NONZERO)=@_;
	
#	print Dumper($IN_FILE);
	
	my $out_file_name = "variance-for-sync-sliding-window.params";
	my ($inName, $path)= fileparse($IN_FILE);
	$out_file_name = $path.$out_file_name;
	
	open fileHandle, ">>", $out_file_name; 
	print "Info about runned parameters is stored in file $out_file_name and next info will be appended to the file if the script will run again.\n";
	print fileHandle "script $RealBin/$RealScript runed at ", POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime), " with the following parameter setting:\n";
	print fileHandle "input sync-variance file: $IN_FILE\n";
	print fileHandle "output file: $out_file_name\n";
	print fileHandle "window-size: $WINDOW_SIZE\n";
	print fileHandle "step-size: $STEP_SIZE\n";
	print fileHandle "min-cov-fraction: $MIN_COV_FRACTION\n";
	print fileHandle "min-length-fraction: $MIN_LENGTH_FRACTION\n";	
	print fileHandle "all-pop-nonzero: $ALL_POP_NONZERO\n\n";
	close fileHandle;
}



=head1 NAME

variance-for-sync-sliding-window.pl

=head1 SYNOPSIS

variance-for-sync-sliding-window.pl --sync-variance INPUT_FILE --output OUTPUT_FILE --window-size 1000 --step-size 1000 --min-cov-fraction 0.6 --min-length-fraction 1 --all-pop-nonzero --print-zero-lines

=head1 DEFAULT SETTINGS

window-size = 1000;
step-size = 1000;
min-cov-fraction = 0.6;
min-length-fraction = 1;
all-pop-nonzero = 0;
print-zero-lines = 0; 
