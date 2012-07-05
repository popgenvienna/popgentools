#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

use Getopt::Long;
use Pod::Usage;

my $HELP;

my $prettybaseFile;
my $syncFile;

GetOptions(
	"sync=s"=>\$syncFile,
	"prettybase-summary=s"=>\$prettybaseFile,
	'help'=> \$HELP,		
) or die "GetOptions";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;


my %insertions=();

my $insLength = 0;

open prettybaseFileHandle, "<", $prettybaseFile;
while (my $line=<prettybaseFileHandle>){
	next if $line =~ m/^#/;
	chomp($line);
	my ($position, $type, $sig, $length) = split "\t", $line;
	
	$insertions{$position}{length} = $length;
	$insertions{$position}{sig} = $sig;
	$insertions{$position}{type} = $type;
}
close prettybaseFileHandle;

open syncFileHandle, "<", $syncFile;
my $line;

$line = <syncFileHandle>;
foreach my $position (sort {$a<=>$b} keys %insertions){

	chomp($line);
	my ($chr,$pos,@parts) = split "\t", $line;

	while ($pos < $position){
		if (!eof(syncFileHandle)){ $line = <syncFileHandle>;}else{last;}
		chomp($line);		
		($chr, $pos, @parts) = split "\t", $line;
	}
		
	print "#".$position."\t".$insertions{$position}{length}."\t".$insertions{$position}{type}."\t".$insertions{$position}{sig}."\n";	
		
	while (($pos >= $position) and ($pos <= $position + $insertions{$position}{length} - 1)){

		print $chr."\t".$pos."\t".join("\t", @parts)."\n";						
		if (!eof(syncFileHandle)){ $line = <syncFileHandle>;}else{last;}
		chomp($line);
		($chr, $pos, @parts) = split "\t", $line;

	}	
	
}

close syncFileHandle;

=pod

=head1 NAME

select-prettybase-inserts-from-sync.pl

=head1 SYNOPSIS

perl select-prettybase-inserts-from-sync.pl --sync fileName --prettybase-summary fileName 

=head1 DESCRIPTION

The script prints those positions of synchronized pileup file that are in insertions of prettybase-summary file. For each insertion/SNP a commented 
line with prettybase-summary is printed and the line is followed by all records from synchronized pileup file that fall into the insertion/SNP region.
The used synchronized pileup file should contain data only for one chromosome.

=head1 OPTIONS

=over 4

=item --sync

A synchronized pileup file. For details about the file format see script synchronize-pileup.pl in popoolation. The synchronized pileup 
file should contain only data for one chromosome.

=item --prettybase-summary

A tab delimited file that contains following columns: position, 
type of feature ("INS" for insertion or "SNP"), significance ("SIG" or "NON"), length. You can create a prettybase-summary file from 
a prettybase file using script PopGenTools/misc/tina/summarize-prettybase.pl, see its help page for more details.
	
=item --help

Prints the help page.

=head1 EXAMPLE

=over 4

=item INPUT example.sync

  2L	1	A	180:10:0:0:0:0
  2L	2	T	0:150:0:0:0:0
  2L	5	G	0:0:15:97:0:0
  2L	18	T	0:139:0:0:0:0
  2L	21	G	0:0:6:93:0:0

=item INPUT example.summary

  1	INS	NON	6

=item Used command

perl select-prettybase-inserts-from-sync.pl --sync example.sync --prettybase-summary example.summary
 
=item OUTPUT

  #1	6	INS	NON
  2L	1	A	180:10:0:0:0:0
  2L	2	T	0:150:0:0:0:0
  2L	5	G	0:0:15:97:0:0
