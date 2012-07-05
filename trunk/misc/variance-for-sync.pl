#!/usr/bin/perl

use strict;
use warnings;

use FindBin qw($Bin $Script);
use lib "$Bin/../Modules";

use Data::Dumper;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX;


use TinaPerl;
 
my $help = 0;
my $test="";

my $SYNC_PILEUP_FILE="";
my $PI_FILE;
my $THETA_FILE;
my $D_FILE;


my $POOL_SIZE = 0;



###default settings:

my $MIN_COUNT = 2;

my $MIN_COV = 4;

my $MAX_COV = 100000;

GetOptions(
	"sync-pileup=s"=>\$SYNC_PILEUP_FILE,
	"out-pi=s"=>\$PI_FILE,
	"out-theta=s"=>\$THETA_FILE,
	"out-D=s"=>\$D_FILE,
	"min-count=i"=>\$MIN_COUNT, 
	"min-cov=i"=>\$MIN_COV,
	"max-cov=i"=>\$MAX_COV,
	"pool-size=i"=>\$POOL_SIZE,
	"test=s"=>\$test,
  	"help"=>\$help
) or die "Invalid arguments, use 'perl variance-for-sync.pl --help'.";

#
# testing:
#

#$POOL_SIZE=100;
#my $STR;
#$STR = 
#"2L\t5002\tG\t0:6:0:32:0:0\t0:4:0:22:0:0\t0:2:0:36:0:0\t0:4:0:23:0:0\n".
#"2L\t5465\tC\t4:0:46:0:0:0\t8:0:66:0:0:0\t0:4:67:0:0:0\t0:4:53:0:0:0\n".
#"2L\t5495\tA\t51:1:6:0:0:0\t49:4:0:0:0:0\t53:6:0:0:0:0\t6:42:0:0:0:0\n".
#"2L\t5762\tT\t0:38:16:0:0:0\t0:19:16:0:0:0\t0:32:23:0:0:0\t0:4:25:0:0:0\n".
#"2L\t5776\tC\t0:0:57:2:0:0\t0:0:33:0:0:0\t0:0:43:0:0:0\t0:0:19:0:0:0\n";

#$SYNC_PILEUP_FILE=\$STR;
#$THETA_FILE="/Users/tinavisnovska/data/Pablo/testTheta";

#
# test end
#

if ($test eq "prove"){
	my $testOut;
	$testOut = `prove $Bin/../Modules/Test/varianceForSync.t`;
	print $testOut, "\n";
}elsif ($test eq "perl"){
	my $testOut;
	$testOut = `perl $Bin/../Modules/Test/varianceForSync.t`;	
	print $testOut, "\n";

}else{

	pod2usage(-verbose=>3) if $help;

	pod2usage(-msg=>"Could not find input sync pileup file.",-verbose=>1) unless defined($SYNC_PILEUP_FILE);

	if (!(defined($PI_FILE))and(!defined($THETA_FILE))and!(defined($D_FILE))){pod2usage(-msg=>"Could not find any output file (specify at least one).",-verbose=>1);}

	if ($POOL_SIZE == 0){pod2usage(-msg=>"Do not specified required option --pool-size.",-verbose=>1)};


	my $fileHandlePi;
	my $fileHandleTheta;
	my $fileHandleD;

	my $ptrOut = [];
	if (defined($PI_FILE)){
		open $fileHandlePi, ">>", $PI_FILE;
		push @{$ptrOut}, {filehandle=>$fileHandlePi,measure=>"pi",path=>$PI_FILE};
	}

	if (defined($THETA_FILE)){
		open $fileHandleTheta, ">>", $THETA_FILE;
		push @{$ptrOut}, {filehandle=>$fileHandleTheta,measure=>"theta",path=>$THETA_FILE};
	}

	if (defined($D_FILE)){
		open $fileHandleD, ">>", $D_FILE;
		push @{$ptrOut}, {filehandle=>$fileHandleD,measure=>"D",path=>$D_FILE};
	}
	
#	print Dumper($SYNC_PILEUP_FILE);

	print_input_params_append($SYNC_PILEUP_FILE, $ptrOut, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV);	
	synchronized_measures_low_memory_requirement($SYNC_PILEUP_FILE, $ptrOut, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV);
	
	
	foreach my $ptrFileInfo (@{$ptrOut}){
		close $ptrFileInfo->{filehandle};	
	}

	

}



sub print_input_params_append{
	my ($IN_FILE, $ptrSettings, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV)=@_;
	
#	print Dumper($IN_FILE);
	
	my $out_file_name = "variance-for-sync.params";
	my ($inName, $path)= fileparse($IN_FILE);
	$out_file_name = $path.$out_file_name;
	
	open fileHandle, ">>", $out_file_name; 
	print "Info about runned parameters is stored in file $out_file_name and next info will be appended to the file if the script will run again.\n";
	print fileHandle "script $Bin/$Script runed at ", POSIX::strftime("%m/%d/%Y %H:%M:%S\n", localtime), " with the following parameter setting:\n";
	print fileHandle "input sync-pileup file: $IN_FILE\n";

	foreach my $ptrSetting (@$ptrSettings){
		print fileHandle "output file: ", $ptrSetting->{path}, "\n";			
	}

	print fileHandle "pool-size: $POOL_SIZE,\n";
	print fileHandle "min-count: $MIN_COUNT\n";
	print fileHandle "min-cov: $MIN_COV\n";
	print fileHandle "max-cov: $MAX_COV\n\n";
	close fileHandle;
}






=head1 NAME

variance-for-sync.pl

=head1 SYNOPSIS

perl variance-for-sync.pl --sync-pileup syncFileName --out-pi outPiFileName --pool-size 100

=head1 OPTIONS

The required options:

=over 10



=item B<--sync-pileup> - An input synchronized pileup file.


=item B<--pool-size> - Size of your pool.


=item B<--out-pi> 

=item B<--out-theta>

=item B<--out-D>

Output files for different measures, at least one of the files should be specified. 



=back

The options that have default values alredy set and if the setting is fine for you there is no need to set them again:

=over 10

=item B<--min-count> - The minimum count of nucleotides of the same type that are needed for SNP identification. 



default = 2



=item B<--min-cov> - The minimum coverage at a SNP site.



default = 4



=item B<--max-cov> - The maximum coverage at a SNP site.


default = 100000


=back

The options that leads to unit testing or help page: 



=over 9



=item B<--test> - Run the unit tests for this script.



=over



=item  --test prove 



Runs prove unit tests, suggested. 



=item --test perl



Runs more detailed perl unit test.



=back



=item B<--help> - Display help for this script.



=back 