#!/usr.bin/perl


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Pod::Usage;


my $FILE_NAME;
my @POPULATIONS_FIRST; 
my @POPULATIONS_SECOND; 
my @ALL_POPS;
my $HELP;

GetOptions(
	"input=s"=>\$FILE_NAME,
	"populations=s"=>\@ALL_POPS,
	"help"=>\$HELP,
)or die "GetOptions";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|EXAMPLE|OPTIONS"}) if $HELP;


@ALL_POPS = (split(/,/, join(',', @ALL_POPS)));

foreach my $pair (@ALL_POPS){
	my ($first, $second) = split "-", $pair;
	push @POPULATIONS_FIRST, $first;
	push @POPULATIONS_SECOND, $second;	
}

open fileHandle, "<", $FILE_NAME;

while (my $line =<fileHandle>){

	next if $line=~m/^#/;
	chomp($line);
	
	print $line;
	
	my ($chr, $pos, $ref, @info) = split "\t", $line;
	
	my $i=0;
	my $set_alleles = 1;
	my $A_index;	my $a_index;
	
	my %SumCounts=();		
	my %counts1=();	
	my $count_A;
	my $count_a;
	
	my $count1_A;
	my $count1_a;
	
	while ($i<scalar(@POPULATIONS_FIRST)){
		
		my $pop1 = $info[$POPULATIONS_FIRST[$i]-4];
		my $pop2 = $info[$POPULATIONS_SECOND[$i]-4];

		my @parts1 = split ":", $pop1;
		my @parts2 = split ":", $pop2;
				
		$SumCounts{0}=$parts1[0]+$parts2[0];
		$SumCounts{1}=$parts1[1]+$parts2[1];
		$SumCounts{2}=$parts1[2]+$parts2[2];
		$SumCounts{3}=$parts1[3]+$parts2[3];
		
		if ($set_alleles){

			$counts1{0}=$parts1[0];
			$counts1{1}=$parts1[1];
			$counts1{2}=$parts1[2];
			$counts1{3}=$parts1[3];
		
			my $ptrSortedByValue = [];
			foreach my $nucleotide (sort {$SumCounts{$b}<=>$SumCounts{$a}} (keys(%SumCounts))){
				my %hash = (val=> $SumCounts{$nucleotide}, nucleotide => $nucleotide);
				push @{$ptrSortedByValue}, \%hash;
			}

			$A_index = $ptrSortedByValue->[0]{nucleotide};
			$a_index = $ptrSortedByValue->[1]{nucleotide};
			
			$count_A = $ptrSortedByValue->[0]{val};
			$count_a = $ptrSortedByValue->[1]{val};		
				
			$count1_A = $parts1[$A_index];
			$count1_a = $parts1[$a_index];
			
			$set_alleles = 0;		
		}else{
			$count_A = $SumCounts{$A_index};
			$count_a = $SumCounts{$a_index};
				
			$count1_A = $parts1[$A_index];
			$count1_a = $parts1[$a_index];
		}
		
		my $odds_ratio;
			
			$odds_ratio	= ($count1_A + 0.5)*(($count_a-$count1_a) + 0.5) / ( ($count1_a+0.5) * (($count_A-$count1_A)+0.5));
			#print "\t".$count1_A."\t".$count1_a."\t".($count_A-$count1_A)."\t".($count_a-$count1_a)."\t".$odds_ratio;
			print "\t".$odds_ratio;	
		$i+=1;	
	}
	print "\n";
}

close fileHandle;


sub max{
	my $m=$_[0];
	for (my $i=0; $i< scalar(@_); $i++){
		next unless ($_[$i]>$m);
		$m = $_[$i];		
	} 
	return $m;
}


=pod

=head1 NAME 

calculate-odds-ratios.pl

=head1 SYNOPSIS

perl calculate-odds-ratios.pl --input syncFileName --populations 4-6,5-7 > outputFileName

=head1 DESCRIPTION 

The script prints input file with some additional columns, 
one for each pair of populations defined by parameter --populations. 
In the additional columns odds ratios are printed for the specified pairs of columns. 
Odds ratios are corrected for continuity. 

=head1 OPTIONS

=over 4

=item --input

A synchronized pileup file.

=item --populations
	
Pairs of populations for which the odds ratio will be calculated.	
Each pair of populations has to be separated by a "," and the two
populations for which the odds ratio will be calculated by a "-". For example when
user provides 4-6,5-7 the script will calculate odds ratio for populations in fourth and sixth column 
and  odds ratio for populations in fifth and seventh column. The columns are numbered from 1. 
Usually, the first population is in column 4. 

=item --help, -h

Prints the help page.

=back

=head1 EXAMPLE

=over 4

=item INPUT file example.sync
	
	2L	1	A	0:20:0:89:0:0	0:15:0:95:0:0	0:110:0:25:0:0	0:95:0:30:0:0
	2L	3	C	0:20:0:89:0:0	0:15:0:95:0:0	0:25:0:113:0:0	0:30:0:111:0:0

=item Used command

perl calculate-odds-ratios.pl --input example.sync --populations 4-6,5-7

=item OUTPUT
	2L	1	A	0:20:0:89:0:0	0:15:0:95:0:0	0:110:0:25:0:0	0:95:0:30:0:0	0.0528577567683713	0.0518352018859132
	2L	3	C	0:20:0:89:0:0	0:15:0:95:0:0	0:25:0:113:0:0	0:30:0:111:0:0	0.98087461050822	1.68537537971937

=back

=cut
