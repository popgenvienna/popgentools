#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $SYNC_FILE_NAME="";
my @COLUMNS=();
my $HELP;

GetOptions(
	"input=s"=>\$SYNC_FILE_NAME,
	"populations=s"=> \@COLUMNS,
	'help'=> \$HELP,
);

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

if (!defined($SYNC_FILE_NAME) or (scalar @COLUMNS == 0)){
  die "Invalid arguments, use 'perl $0 --help'.";
}

open fileHandle, "<", $SYNC_FILE_NAME or die "Could not open file $SYNC_FILE_NAME";

@COLUMNS = (split(/,/, join(',', @COLUMNS)));

while (my $line = <fileHandle>){
	if ($line =~ m/^#/){
	  print $line;
	}else{
	  chomp($line);
		my @parts = split "\t", $line;		
		my @coverages=();
  
		for (my $i=0; $i<= (scalar(@COLUMNS) -1); $i++ ){

		  my ($a,$t,$c,$g) = split ":", $parts[$COLUMNS[$i]-1];			
		  push @coverages, $a+$t+$c+$g;
    
		}
	
	  my $min = min(@coverages);
	
	  print $line."\t".$min."\n";
	}
}

sub min{
	my $m=$_[0];
	for (my $i=0; $i< scalar(@_); $i++){
		next unless ($_[$i]<$m);
		$m = $_[$i];		
	} 
	return $m;
}

=pod

=head1 NAME 

add-min-coverage.pl

=head1 SYNOPSIS

perl add-min-coverage.pl --input fileName --populations 4,5

=head1 DESCRIPTION

The script takes an input synchronized pileup file and calculates minimal coverage over population columns 
specified by option --populations.

=head1 OPTIONS

=over 4

=item --input

A synchronized pileup file.

=item --populations

A comma separated list of column numbers for columns that contains populational data over which the minimal 
coverage should be calculated. A column with a chromosome name has a number 1, first populational column 
has a number 4.

=item --help, -h

Prints the help page.

=back

=head1 EXAMPLE

=over 4

=item INPUT example.sync

  2L	1	G	0:1:15:90:0:0	0:0:25:100:0:0	0:0:80:20:0:0	0:0:95:15:0:0
  2L	6	C	0:1:15:90:0:0	0:0:25:100:0:0	0:0:80:20:0:0	0:0:95:15:0:0
  2L	10	C	0:0:80:1:0:0	0:0:115:0:0:0	0:1:80:3:0:0	0:1:95:0:0:0  

=item Used command

perl add-min-coverage.pl --input example.sync --populations 4,5,6,7

=item OUTPUT

  2L	1	G	0:1:15:90:0:0	0:0:25:100:0:0	0:0:80:20:0:0	0:0:95:15:0:0	100
  2L	6	C	0:1:15:90:0:0	0:0:25:100:0:0	0:0:80:20:0:0	0:0:95:15:0:0	100
  2L	10	C	0:0:80:1:0:0	0:0:115:0:0:0	0:1:80:3:0:0	0:1:95:0:0:0	81

