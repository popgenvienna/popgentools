#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use Synchronized;



# --pool-size 500 --min-count 2 --min-coverage 4 --window-size 1000 --step-size 1000 --input test/snp.merge --output test/test.fst
my $input;
my $output;
my $help=0;
my $test=0;
my $mincount=0;
my $mincoverage=0;
my $maxcoverage=0;
my $chromosomes="";


GetOptions(
    "input=s"	        =>\$input,
    "output=s"          =>\$output,
    "min-count=i"       =>\$mincount,
    "min-coverage=i"    =>\$mincoverage,
    "max-coverage=i"    =>\$maxcoverage,
    "chromosomes=s"     =>\$chromosomes,
    "test"              =>\$test,
    "help"	        =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);
pod2usage(-verbose=>2) if $help;
pod2usage() unless -e $input;
pod2usage() unless $output;
pod2usage() unless $mincoverage;
pod2usage() unless $mincount;
pod2usage() unless $maxcoverage;

my @temp=split /\s+/,$chromosomes;
push @temp, $chromosomes unless @temp;
my %chrh=map {($_,1)} @temp;

my $sp=get_sumsnp_synparser($mincount,$mincoverage,$maxcoverage);

open my $ifh,"<", $input or die "Could not open input file";
open my $ofh,">",$output or die "Could not open output file";

while(my $line=<$ifh>)
{
    my $p=$sp->($line);
    my $chr=$p->{chr};
    next unless(exists($chrh{$chr}));
    next unless $p->{ispuresnp};
    print $ofh $line;
}
exit;

=head1 NAME

filter-synchronized.pl - Filter a synchronized file by some  

=head1 SYNOPSIS

 filter-synchronized.pl --input populations.sync --output filtered.sync --min-count 2 --min-coverage 4  --max-coverage 1000000 --chromosomes "2L 2R X 3R 3L 4"
 
=head1 OPTIONS

=over 4

=item B<--input>

The input file. Has to be synchronized pileup file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--chromosomes>

List of chromosomes which should be used; eg: "2L 2R X 3L 3R 4"

=item B<--min-count>

the minimum count of the minor allele. used for SNP identification. SNPs will be identified considering all populations simultanously. 

=item B<--min-coverage>

the minimum coverage; used for SNP identification, the coverage in ALL populations has to be higher or equal to this threshold, otherwise no SNP will be called. 

=item B<--max-coverage>

the maximum coverage; used for SNP identification, the coverage in ALL populations has to be lower or equal to this threshold, otherwise no SNP will be called. 

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=cut





