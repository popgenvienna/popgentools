#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $HELP;

my $positions;
my $input;
my $filter=2;

GetOptions(
	"positions=s"=>\$positions,
	"input=s"=>\$input,	
	"filter=i"=>\$filter,
	'help'=> \$HELP,
) or die "GetOptions";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

if (!defined($input) or !defined($positions)){die "Options --positions and --input should be specified."}

if (($filter != 1) and ($filter != 2)){die "Wrong --filter parameter value."}

my %filterPos;

open POSITIONS_FILE_HANDLE, "<", $positions or die "Could not open a file with positions.";

while (my $line = <POSITIONS_FILE_HANDLE>){
 	next if $line =~ m/^#/;
 	chomp($line);
 	my ($chromosome, $position) = split "\t", $line;

 	$filterPos{$chromosome}{$position}=1;
}
close POSITIONS_FILE_HANDLE;


open INPUT_FILE_HANDLE, "<", $input or die "Could not open a pileup-sync-pileup file.";
while (my $line = <INPUT_FILE_HANDLE>){
	next if $line=~m/^#/;
	chomp($line);
	my ($chromosome, $position) = split "\t", $line;
	
	if ($filter == 1){
		#keep
		if (exists($filterPos{$chromosome}{$position}) and ($filterPos{$chromosome}{$position} == 1)){
			print $line."\n";
		}	
	}elsif( $filter == 2){
		#filter out
		if (!(exists($filterPos{$chromosome}{$position}) and ($filterPos{$chromosome}{$position} == 1))){
			print $line."\n";
		}			
	};
	
}
close INPUT_FILE_HANDLE; 


1;

__END__

=pod

=head1 NAME

filter-positions-pileup.pl

=head1 SYNOPSIS

perl filter-positions-pileup.pl --input inputFileName --positions posFileName --filter 1 

=head1 DESCRIPTION

The script takes pileup or sync-pileup file and filter out (or keep) the positions specified in a file positions.

=head1 OPTIONS

=over 4

=item --input

An input file, pileup or sync pileup or output of CMH-test, any tab delimited file that contain chromosome and position in the first two columns.

=item --positions

A tab delimited file used for filtering, consists of two columns, the first is chromosome and the second is position on chromosome.

=item --filter 

A value of the parameter should be 1 or 2. If the value is 1 the output contains those lines from --input file which positions are specified in --positions. If the value is 2 all the specified positions are filtered out from --input file. (default value 2) 

=item --help, -h

Prints the help page.

=back

=head1 EXAMPLE

=over 4

=item INPUT file example.sync

	2L	1	0:0:15:90:0:0
	2L	5	100:0:0:0:0:0
	3L	18	1:25:70:0:0:0
	3L	25	1:0:80:15:0:0
	
=item INPUT file positions.list

	2L	5
	3L	18
	
=item Used command 1

perl filter-positions-pileup.pl --input example.sync --positions positions.list --filter 1

=item OUTPUT 1

	2L	5	100:0:0:0:0:0
	3L	18	1:25:70:0:0:0

=item Used command 2

perl filter-positions-pileup.pl --input example.sync --positions positions.list --filter 2

=item OUTPUT 2

	2L	1	0:0:15:90:0:0
	3L	25	1:0:80:15:0:0

=back
	
=cut
