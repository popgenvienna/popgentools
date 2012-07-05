#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $prettybaseSummarized;
my $bam;
my $samtools;

my $HELP;
  
GetOptions(
	"input=s"=>\$prettybaseSummarized,
	"bam=s"=>\$bam,
	"samtools=s"=>\$samtools,
	'help'=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

my %insertions=();


print "#chromosome\tposition\ttype\tlength\tstart_coverage\tend_coverage\n";

open prettybaseSummarizedFileHandle, "<", $prettybaseSummarized;
while (my $line = <prettybaseSummarizedFileHandle>){
	next if $line =~ m/^#/;	
	chomp($line);
	my ($chromosome, $position, $type, $length) = split "\t", $line;
	if ($type eq "INS"){
	#print coverage of breakpoints	
		my $start = `$samtools view -c $bam $chromosome:$position-$position`;
		chomp($start);
		my $endpos = $position + $length - 1;
		my $end = `$samtools view -c $bam $chromosome:$endpos-$endpos`; 
		chomp($end);
	#	print "start_cov\t".$start."\tend_cov\t".$end."\n";
		print $chromosome."\t".$position."\t".$type."\t".$length."\t".$start."\t".$end."\n";
	}
}

close prettybaseSummarizedFileHandle;

=pod

=head1 NAME 

get-breakpoints-coverage.pl

=head1 SYNOPSIS

perl get-breakpoints-coverage.pl --input prettybaseSummaryFile --bam bamFile --samtoolsBin samtoolsBinary

=head1 DESCRIPTION

The script takes an input file in prettybase-summarized format and a bam file and for each region defined in prettybase-summarized format 
the script returns a coverage of start and end position. 

Output is a tab separated table with the following columns: chromosome, position, type (INS or SNP, the same as in the input 
prettybase-summarized file), length (as in the input), start position coverage, end position coverage. 

=head1 OPTIONS

=over 4

=item --input 

A file in prettybase-summarized file format (for details see script summarize-prettybase.pl).

=item --bam

A bam file that contains information about mapping to the same region as from witch the prettybase-sumarized file is.

=item --samtools

A path to a binary of samtools, version that allows you to specify samtools view -c option.

=item --help, -h

Prints the help page. 
