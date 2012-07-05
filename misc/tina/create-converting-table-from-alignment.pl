#!/usr/bin/perl

use warnings;
use strict; 
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $HELP;

GetOptions(
	'help'=> \$HELP,
) or die "GetOptions";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|EXAMPLE|OPTIONS"}) if $HELP;

#my $alignmentFileName = shift;

my $usage = "perl $0";
# print convert table from 	

#open fileHandle, "<", $alignmentFileName or die $usage;

my $ptrHeaders=[];
my $ptrSequences=[];

while (my $line = <STDIN>){
	
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

my $length = length($ptrSequences->[0]);

print "#alignment length: ".$length."\n";
print "#headers: ".$ptrHeaders->[0].", ".$ptrHeaders->[1]."\n";
my @gaps = (0,0);

print "#position_in_reference_of_bab_sequence_0=1029080\tbase_in_ref_bab_sequence\tposition_in_bickels_bab_sequence_0based\tbase_in_bickels_bab_sequence\n";

my $babPrevPosition = 1;
my $refPrevPosition = 1;

for (my $i=0; $i< $length; $i++){
	my $printout=1;	

	for (my $j=0; $j<=1; $j++){
		if (substr($ptrSequences->[$j], $i, 1) eq "-"){
			$gaps[$j]+=1;	
			$printout=0;
		}
	}

#	print Dumper($printout);	
	if ($printout){
		my $bickelPosition = $i - $gaps[0];
		my $bickelLetter = substr($ptrSequences->[0], $i, 1);
	
		my $refPosition = $i - $gaps[1];
		my $refLetter = substr($ptrSequences->[1], $i,1);
		
#		my $refPrint = $refPosition+1029080-1;

		print $refPosition."\t".$refLetter."\t".$bickelPosition."\t".$bickelLetter."\n";
	}
}
	
1;

__END__

=pod

=head1 NAME

create-converting-table-from-alignment.pl

=head1 SYNOPSIS

perl create-converting-table-from-alignment.pl < fastaFileName 

=head1 DESCRIPTION

The script requires an input fasta file where a pairwise alignment is stored. The script creates a converting table between base numberings of sequences aligned in the input file. Output is tab delimited file with numbers of positions that correspond to each other in one line. 

Output consists of four columns: position_number_in_the_second_sequence, base_on_the_position_of_the_second_sequence, position_number_in_the_first_sequence, base_on_the_position_of_the_first_sequence, positions are 0 based (the first position has number 0).


=head1 OPTIONS

=over 4 

=item --help, -h

Prints the help page.

=back

=head1 EXAMPLE

=over 4

=item INPUT file example.fasta

	>seq1
	ATCG--CAGTAC
	>seq2
	A-CGAAC--TAT


=item Used command

perl < example.fasta 

=item OUTPUT 

	0	A	0	A
	1	C	2	C
	2	G	3	G
	5	C	4	C
	6	T	7	T
	7	A	8	A
	8	T	9	C

=back

=cut
