#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

my $HELP;

my $FASTA_FILE;
my $DIVERGENCE_RATIO;
my $SIMILARITY;

GetOptions(
	"input=s"=>\$FASTA_FILE,
	"similarity=i"=>\$SIMILARITY,
	'help'=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

$DIVERGENCE_RATIO=1-$SIMILARITY;

open fileHandle, "<", $FASTA_FILE;
srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);	

while (my $line = <fileHandle>){
	chomp($line);
	next if $line=~m/^#/;
	if ($line =~m/^>/){
		print $line, "\n";	
	}else{
		#srand(10);	

		my $tmp="";
		
		my $c=substr($line, 0, 1);
		$line = substr($line, 1);
		
		while(length($line)>=1){
			my $r = rand();
			
			$c = randomly_diverge($c);
			print $c;
				
			$c=substr($line, 0,1);
			$line = substr($line, 1);	
		}

		$c= randomly_diverge($c);
		print $c;
		
		print "\n";
	}
}
close fileHandle;
		
		
sub randomly_diverge{
	my ($in_c)=@_;	
	my $out_c = $in_c;	
	
	my $r = rand();
			
	if ($r <= $DIVERGENCE_RATIO){
	#change $c randomly to another base
	my $s = rand();			
	if($out_c eq "A"){
		if($s <= 1/3){$out_c = "T"}
		elsif((1/3<=$s)and($s<=2/3)){$out_c="C"}
		else{$out_c="G"}
	}elsif($out_c eq "T"){
		if($s <= 1/3){$out_c = "A"}
		elsif((1/3<=$s)and($s<=2/3)){$out_c="C"}
		else{$out_c="G"}
	}elsif($out_c eq "C"){
		if($s <= 1/3){$out_c = "T"}
		elsif((1/3<=$s)and($s<=2/3)){$out_c="A"}
		else{$out_c="G"}
	}elsif($out_c eq "G"){
		if($s <= 1/3){$out_c = "T"}
		elsif((1/3<=$s)and($s<=2/3)){$out_c="A"}
		else{$out_c="C"}	
	}							
	}
	return $out_c;
}		
		

=pod

=head1 NAME 

create-similar-sequence.pl

=head1 SYNOPSIS

perl create-similar-sequence.pl --input fastaFile --similarity 0.6 

=head1 DESCRIPTION

The script takes an input fasta file and for each sequence of the file creates a random sequence which similarity to the original sequence is 
defined by --similarity option.

=head1 OPTIONS

=over 10

=item --input 

An input fasta file.

=item --similarity

A number between 0 and 1. If the value is 1, the generated sequence will be the same as input one, if the value is 0, 
the generated sequence will be completly different.

=head1 EXAMPLE

=over 4

=item INPUT example.fasta

    >seq1
    ATATCTGCGATC
    >seq2 
    CGCTCGATCA

=item Used command

perl create-similar-sequence.pl --input example.fasta --similarity 1

=item OUTPUT

    >seq1
    ATATCTGCGATC
    >seq2 
    CGCTCGATCA

=back

=cut