#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

use Getopt::Long;
use Pod::Usage;

my $HELP;

GetOptions(
	'help'=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP; 

my $line = <STDIN>;
 
my $a = 0;
my $t = 0;
my $c = 0;
my $g = 0;
 
while ($line = <STDIN>){
	chomp($line);
	next if $line =~ m/^>/;	
	
	my @F = split "", $line;
	
	foreach my $f (@F){
		if (($f eq 'a')or($f eq 'A')) {
			$a+=1;	
		}elsif(($f eq 't')or($f eq 'T')){
			$t+=1;	
		}elsif(($f eq 'c')or($f eq 'C')){
			$c+=1;
		}elsif(($f eq 'g')or($f eq 'G')){
			$g+=1;
		}else{
			print STDERR "problem\n";	
		}
	}
}

print "count a: $a\n".
	  "count t: $t\n".
	  "count c: $c\n".
	  "count g: $g\n".
	  "sum : ", $a+$t+$c+$g, "\n";  	

1;

__END__

=pod

=head1 NAME

count-nucleotides.pl

=head1 SYNOPSIS

perl count-nucleotides.pl < inputFastaFile

=head1 DESCRIPTION

The script expects a file with one fasta sequence as an input, counts a numbers of 'a', 't', 'c', and 'g' in the input file and prints out the counts. 

=head1 EXAMPLE

=over 4 

=item INPUT example.fasta

	>seq1
	ATCTGATCGG

=item Used command

perl count-nucleotides.pl < example.fasta

=item OUTPUT

	count a: 2
	count t: 3
	count c: 2
	count g: 3
	sum : 10

