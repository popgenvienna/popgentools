#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $HELP;

my $positionOfSequence=1;
my $listOfPositions;
my $fasta;
my $length=10;


GetOptions(
	"position-of-sequence=i" => \$positionOfSequence,
	"list-of-positions=s" => \$listOfPositions,
	"fasta-pairwise-alignment=s" => \$fasta,
	"flanking-region-length=i" => \$length,
	'help'=> \$HELP,
) or die "GetOptions";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

if (!defined($fasta)){die "Option --fasta-pairwise-alignment not defined."}

if (!defined($listOfPositions)){die "Option --list-of-positions not defined."}


if (($positionOfSequence>1) or ($positionOfSequence<0)){die "Option --position-of-sequence should be 0 or 1."}

my %positions = ();

#read the positions that should be printed with their flanking regions
open listFileHandle, "<", $listOfPositions or die "Could not open a list of positions $listOfPositions";
while (my $line = <listFileHandle>){
	next if $line =~ m/^#/;
	chomp($line);
	$positions{$line} = 1;
}
close listFileHandle;
 

#read fasta alignment
open fileHandle, "<", $fasta or die "Could not open a fasta file $fasta";

my $ptrHeaders=[];
my $ptrSequences=[];
my $ptrGaps = [];

while (my $line = <fileHandle>){
	
	next if $line =~ m/^#/; 
	chomp($line);
	
	if ($line =~ m/^>/){
		push @$ptrHeaders, substr($line,1);
		push @$ptrSequences, ""; 
	}else{
		my $tmpSeq = pop(@$ptrSequences);
		$tmpSeq = $tmpSeq.$line;
		push @$ptrSequences, $tmpSeq;
	}
}
close fileHandle;

my $l = length($ptrSequences->[0]);

for (my $i = 0; $i < $l; $i++){
	
	if (substr($ptrSequences->[0], $i, 1) eq "-"){
		if ($i>=1){
			$ptrGaps->[0][$i]=$ptrGaps->[0][$i-1]+1; 
		}else{
			$ptrGaps->[0][$i]=-1;	
		}
	}else{
		if ($i>=1){
			$ptrGaps->[0][$i]=$ptrGaps->[0][$i-1];	
		}else{
			$ptrGaps->[0][$i]=0;		
		}	
	}	

	if (substr($ptrSequences->[1], $i, 1) eq "-"){
		if ($i>=1){
			$ptrGaps->[1][$i]=$ptrGaps->[1][$i-1]+1;
		}else{
			$ptrGaps->[1][$i]=-1;
		}
	}else{
		if ($i>=1){
			$ptrGaps->[1][$i]=$ptrGaps->[1][$i-1];
		}else{
			$ptrGaps->[1][$i]=0;			
		} 
	}		
}

my $ptrRealPos = [];

for (my $i = 0; $i < $l; $i++){
	$ptrRealPos->[0][$i] = $i - $ptrGaps->[0][$i];
	$ptrRealPos->[1][$i] = $i - $ptrGaps->[1][$i];
}


print "#numberingSelected_from0\tnumberingOther_from0\tbaseSelected\tbaseOther\tsurroundingRegionSelected\tsurroundingRegionOther\n"; 

my $sequenceLength = scalar @{$ptrRealPos->[1-$positionOfSequence]};
#search for positions from a list
foreach my $position (sort {$a<=>$b} keys %positions){
	
	my $p=0;
		
	while (($ptrRealPos->[$positionOfSequence][$p]) < ($position-$length)){
		$p+=1;	
	}

	my $flankingSelected="";
	my $flankingOther="";
	my $otherFirst = $ptrRealPos->[1-$positionOfSequence][$p];

	my $positionSelected=0;
	my $positionOther=0;
	my $baseSelected=0;
	my $baseOther=0;

	my $first = 1;
	my $pos=0;
	my $tmp = $p;
	

	
	while ( ($p<=$sequenceLength-1) and (($ptrRealPos->[$positionOfSequence][$p] <= $position+$length) or ($ptrRealPos->[1-$positionOfSequence][$p] < $otherFirst+2*$length+1))){
		
		my $selected = substr($ptrSequences->[$positionOfSequence], $p, 1);
		my $other = substr($ptrSequences->[1-$positionOfSequence], $p, 1);
		
		if (($first)and($ptrRealPos->[$positionOfSequence][$p] == $position)){

			if ($selected eq "-"){
				$positionSelected = "na";
			}else{
				$positionSelected = $ptrRealPos->[$positionOfSequence][$p];		
			}
			$baseSelected = $selected;
			$flankingSelected = $flankingSelected.uc($selected);
						
			$positionOther = $ptrRealPos->[1-$positionOfSequence][$p];	
			$baseOther = $other;
			$flankingOther = $flankingOther.uc($other);
			
			$first = 0;
			
			$pos = length($flankingSelected);
			
		}else{
			$flankingSelected = $flankingSelected.lc($selected);			
			$flankingOther = $flankingOther.lc($other);			
		}
		$p=$p+1;
	
	}

	$flankingSelected = substr($flankingSelected, $pos-$length-1, 2*$length+1);
	$flankingOther = substr($flankingOther,  $pos-$length-1, 2*$length+1);				
					
	print $positionSelected."\t".$positionOther."\t".$baseSelected."\t".$baseOther."\t".$flankingSelected."\t".$flankingOther."\n";
	

}


=pod

=head1 NAME

flanking-regions.pl

=head1 SYNOPSIS

perl flanking-regions.pl --fasta-pairwise-alignment fileName --list-of-positions listFileName --position-of-sequence 0 --offset 15 --flanking-region-length 4

=head1 DESCRIPTION

The script prints flanking regions for some given positions of a genomic sequence stored in pairwise alignment fasta file.  

For a given position of a given sequence the script prints a specified number of bases before the position and the same number of bases after. The script prints also a corresponding region of the second sequence from the input pairwise alignment fasta file.
A base at the specified position is printed out in upper case, the bases around are printed in lower case.   

Output is a tab separated table with following columns:
	position on selected sequence, 
	position on the other sequence, 
	base of selected sequence, 
	base of the other sequence, 
	flanking refion for selected sequence, 
	flanking region for the other sequence.

=head1 OPTIONS

=over 4

=item --list-of-positions

A path to a file that contains a list of positions of interest. The file contains a position per line. Positions are 0 based (the first base of sequence is on position 0).

=item --fasta-pairwise-alignment

A path to a file that contains a pair of aligned sequences that will be used to identify flanking regions.  

=item --position-of-sequence

A position of a sequence to which a list of positions reffer within a fasta pairwise alignment file (a value can be either 0 or 1).

=item --flanking-region-length

A number of bases that will be printed before and also after the selected position for each flanking region (default value set to 10).

=item --help,-h

Prints the help page.

=back

=head1 EXAMPLE

=over 4

=item INPUT file listOfPositions

	5
	17

=item INPUT file fastaPairwiseAlignment.fasta

	>seq1
	ACAGCTGCTTGGTA-AACCAC
	>seq2
	A-TGCTGCTAGCTAGTACC-C
		
=item Used command

perl flanking-regions.pl --list-of-positions listOfPositions --fasta-pairwise-alignment fastaPairwiseAlignment.fasta --position-of-sequence 0 --flanking-region-length 4 > output.flanking
		
=item OUTPUT file output.flanking:

	5	4	T	T	cagcTgctt	-tgcTgcta
	17	17	C	C	-aacCac	gtacC-c
	
=back	

=cut








