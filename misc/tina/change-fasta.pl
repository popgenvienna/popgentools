#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $HELP;

my $start;
my $end;
my $header;
my $input;
my $op="del";

GetOptions(
	"start=i"=>\$start,
	"end=i"=>\$end,
	"seq-header=s"=>\$header,
	"input=s"=>\$input,
	"operation=s"=>\$op,
	'help'=>\$HELP,
) or die "GetOptions";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

my $usage = "perl $0 header(without >) start end < infile > outfile";


my $selected_seq=0;
my $subseq=""; 

open fileHandle, "<", $input or die "Could not open input file $input.";

while (my $line = <fileHandle>){
	next if $line=~m/^#/;
	chomp($line);
	
	if ($line=~ m/^>/){
		my $line2 = substr($line, 1);
		if ($line2 eq $header){
			$selected_seq = 1;	
		}else{
			$selected_seq = 0;	
			print $line."\n";
		}
	}else{
		if ($selected_seq == 1){
			$subseq=$subseq.$line;	
		}else{
			print $line."\n";	
		}	
	}  
}

close fileHandle;

print ">".$header."\n";

if ($op eq "del"){
	print substr($subseq, 0, $start-1);
	print substr($subseq, $end+1);
}elsif($op eq "sel"){
	print substr($subseq, $start, $end-$start+1);	
}

1;

__END__

=pod

=head1 NAME

change-fasta.pl

=head1 SYNOPSIS

perl change-fasta.pl --operation del/sel --input fileName --seq-header dmel_2L --start 15 --end 158

=head1 DESCRIPTION

The script takes input fasta file and a header of the sequence of interest 
and delete or select a region of the fasta sequence specified by 
--start and --end position parameters.

=head1 OPTIONS

=over 4

=item --operation 

A parameter that specifies which operation takes place. Should be either "del" (deletion) or "sel" (selection). Default value is "del".

=item --input 

An input fasta file.

=item --seq-header 

A header of the sequence from which you want to delete (or select) a region. The header should be without '>' character at the beginning. 

=item --start

The first position of the region that should be deleted (selected) from the fasta sequence. First nucleotide of fasta sequence is on position 0. 

=item --end 

The last position of the region that should be deleted (selected) from the fasta sequence. First nucleotide of fasta sequence is on position 0.

=back

=head1 EXAMPLE 

=over 4

=item INPUT file example.fasta

	>seq1
	attacagatgcgcta
	>seq2
	cgatcgacccgattaga

=item Used command 1

perl change-fasta.pl --operation del --input example.fasta --seq-header seq2 --start 5 --end 10

=item OUTPUT 1

	>seq2
	cgatcattaga

=item Used command 2

perl change-fasta.pl --operation sel --input example.fasta --seq-header seq2 --start 5 --end 10

=item OUTPUT 2

	>seq2
	gacccg
	
=back

=cut