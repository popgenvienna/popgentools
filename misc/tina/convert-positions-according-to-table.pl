#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $INPUT="";
my $CHROMOSOME;
my $CONVERT_TABLE="";
my $DIRECTION;
# 1 for converting from bab_reference to bickels positions
# -1 for converting from bickels to bab_reference positions 
my $HELP;

GetOptions(
	"input=s"=>\$INPUT,
	"convert-table=s"=> \$CONVERT_TABLE,
	"chromosome=s"=> \$CHROMOSOME,
	"direction=i"=> \$DIRECTION,
	'help'=> \$HELP,
) or die "GetOptions";


pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

my %hash=();

open handleConvert, "<", $CONVERT_TABLE or die "Could not open a convert table file: ".$CONVERT_TABLE.".";

if ( ! ($DIRECTION == 1 || $DIRECTION == -1) ){die "Direction should be either 1 or -1."}

if ($DIRECTION == 1){
	#store info in a hash, key - the first column numbering, 	
	while (my $line = <handleConvert>){
		next if $line =~ m/^#/;
		chomp($line);
	
		my ($numbering1, $base1, $numbering2, $base2) = split "\t", $line;
	
		$hash{$numbering1}{numbering2}=$numbering2;
		$hash{$numbering1}{base2}=$base2;
		$hash{$numbering1}{base1}=$base1;		
	}
	
}elsif($DIRECTION == -1){
	#store info in a hash, key - the third column numbering, 	
	while (my $line = <handleConvert>){
		next if $line =~ m/^#/;
		chomp($line);
	
		my ($numbering1, $base1, $numbering2, $base2) = split "\t", $line;

		$hash{$numbering2}{numbering1}=$numbering1;
		$hash{$numbering2}{base1}=$base1;
		$hash{$numbering2}{base2}=$base2;
	}
}
close handleConvert;


open handleInput, "<", $INPUT or die "Could not open an input file ".$INPUT.".";

print "# positions of input file ".$INPUT." converted according to convert table ".$CONVERT_TABLE.", direction ".$DIRECTION."\n";
	
if ($DIRECTION == 1){
	#convert from the first column numbering to the third column numbering	
	while (my $line = <handleInput>){
		if ($line =~ m/^#/){print $line;}
		chomp($line);
		my ($chromosome, $position, $base, @rest) = split "\t", $line;
		
		if ($chromosome eq $CHROMOSOME && defined($hash{$position})){
			print $chromosome."\t".$hash{$position}{numbering2}."\t".$hash{$position}{base2}."\t".join("\t", @rest)."\n";	
		}	
	}	
}elsif($DIRECTION == -1){
	#convert from the third column numbering to the first column numbering		
	while (my $line = <handleInput>){
		if ($line =~ m/^#/){print $line;}
		chomp($line);
		my ($chromosome, $position, $base, @rest) = split "\t", $line;

		if ($chromosome eq $CHROMOSOME && defined($hash{$position})){
			print  $chromosome."\t".$hash{$position}{numbering1}."\t".$hash{$position}{base1}."\t".join("\t", @rest)."\n";
		}	
	}			
}
close handleInput;

__END__

=pod

=head1 NAME

convert-positions-according-to-table.pl

=head1 SYNOPSIS

perl convert-positions-according-to-table.pl --input inputFileName --convert-table convertTableFileName --chromosome 3L --direction 1 

=head1 DESCRIPTION

The script takes an input file and exchanges some positions numbers and nucleotide bases according to a convert table. 

A convert table is a file that contains two different numberings of one region and the numberings 
are based on pairwise alignment of two different sequences from the region (for more details see option "--convert-table"). 
For one convert table, two directions of conversion are possible: from the first numbering to the second one 
(--direction 1) or vice versa (--direction -1). Therefore, the direction of conversion needs to be specified by option "--direction". 

An input file is a tab delimited file that contains chromosome, position, and base in the first three columns respectively. 
The script process the input file line by line. Each time the script finds a line for a position that is defined in a convert 
table, the position and base are exchanged by the corresponding position and base from the convert table.   

=head1 OPTIONS

=over 4

=item --input

A tab delimited file that contains chromosome, position, and base in the first three columns respectively, mandatory

=item --convert-table 

A tab delimited file that contains two different numberings of bases in a region. In the first column, the first numbering position is stored, in the second column, a base that occured on the position. In the third column, the second numbering position is stored and in the fourth column, a base that occured on the position. The data of "--convert-table" file should be only from one chromosome. Data like this can be created by script create-converting-table-from-alignment.pl, see details on the script help page. This option is mandatory for the script.

=item --direction

Set 1 if you want to convert positions from the first numbering to the second one, set -1 if you want to convert positions from second numbering to the first one, mandatory

=item --chromosome

In one run of the script, you can convert only data from one chromosome. In the case when "--input" file contains more than one chromosome the script would translate also positions that are on other chromosomes. To avoid this, you have to specify option "--chromosome". The option is mandatory.

=item --help, -h

Prints the help page.

=back

=head1 EXAMPLE

=over 4

=item INPUT file example.input.sync

	3L	7	T	0:112:13:0:0:0	0:97:5:0:0:0
	3L	9	C	0:20:100:0:0:0	0:16:102:0:0:0

=item CONVERT TABLE file example.convert.table

	0	A	0	A
	2	C	1	C
	3	G	2	G
	4	C	5	C
	7	T	6	T
	8	A	7	A
	9	C	8	T

=item Used command 

perl convert-positions-according-to-table.pl --input example.input.sync --convert-table example.convert.table --chromosome 3L --direction 1 > example.output.sync

=item OUTPUT example.output.sync

	3L	6	T	0:112:13:0:0:0	0:97:5:0:0:0
	3L	8	T	0:20:100:0:0:0	0:16:102:0:0:0

=back

=cut