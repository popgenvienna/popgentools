#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

my $HELP;

GetOptions(
	'help'=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({
  -exitval=>2,
  -verbose=>99, 
  -sections=>"NAME|SYNOPSIS|DESCRIPTION|EXAMPLE",
}) if $HELP;


my $position = 0;

while (my $line = <STDIN>){
	next if $line =~ m/^#/;
	chomp($line);

	if ($line =~ m/^>/){
	  print $line, ", unknown positions: \n"; 
	  $position=0;
	}else{
	  my @F = split "", $line;	
	  foreach my $f (@F){
	    if ($f eq 'N'){
	      print $position."\n";	
	    }
	    $position+=1;
	  }
	}
}

__END__

=pod

=head1 NAME

print-N-positions.pl

=head1 SYNOPSIS

perl print-N-positions.pl < inputFastaFile

=head1 DESCRIPTION

The script prints positions on which unknown (N) nucleotide occured. The script prints these 
positions for each sequence of an input fasta file.

First nucleotide of a fasta sequence is on position 0.

=head1 EXAMPLE

=over 4

=item INPUT example.fasta

  >seq1
  ATANNNTCTGCNGATC
  >seq2 
  CGCTCNNGATCA

=item Used command 
  
perl print-N-positions.pl < example.fasta

=item OUTPUT

  >seq1, unknown positions: 
  3
  4
  5
  11
  >seq2 , unknown positions: 
  5
  6

=back

=cut
