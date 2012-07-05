#!/usr/bin/perl;

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my $prettybase;
my $list="";
my $desiredSum=192;
my $HELP;

GetOptions(
	"input=s"=>\$prettybase,
	"list-of-significant=s"=>\$list,
	"number-of-individuals=i"=>\$desiredSum,
	"help"=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|OPTIONS|DESCRIPTION|EXAMPLE"}) if $HELP;

#my $usage = "perl $0 prettybase listOfSignificant noIndividuals";
#my $prettybaseFile = shift or die $usage;
#my $list = shift or die $usage;    
  
 
my %significant=();

if ($list ne ""){
	open listFileHandle, "<", $list or die "Could not open file ".$list."\n";
	while (my $line = <listFileHandle>){
		next if $line =~ m/^#/;
		chomp($line);
		$line =~ s/^0*//;
		$significant{$line} = 1;
	}
	close listFileHandle;
}

# read data for one position
# process
# print the position data

my $actPosition;
my %actPositionCounts=('A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, 'N'=>0, 'X'=>0, '-'=>0); 
 
print "#position\ttype\tsignificance\tmaxLengthOfInsertion\talleles\tcounts\n";
 
open prettybaseFileHandle, "<", $prettybase or die "Could not open file ".$prettybase;

my $line = <prettybaseFileHandle>;
chomp($line);
my ($pos, $flyID, $base1, $base2)	 = split "\t", $line;
$pos =~ s/^0*//;

$actPosition = $pos;
$actPositionCounts{$base1}+=1;
$actPositionCounts{$base2}+=1;

my $sum;
my $noOfPositions = 0;
my $noOfIndels = 0;
my $type = "";
my $indelMaxLength=0;
my $n;
my $dash;
my $keys;
my $values;
my $indel;


if (length($base1)>=$indelMaxLength){	$indelMaxLength = length($base1);	};
if(length($base2)>=$indelMaxLength){	$indelMaxLength = length($base2);	};




while ($line = <prettybaseFileHandle>){
	chomp($line);
	my ($pos, $flyID, $base1, $base2) = split "\t", $line;
	$pos =~ s/^0*//;

	if ($pos == $actPosition){
		$actPosition = $pos;
		$actPositionCounts{$base1}+=1;
		$actPositionCounts{$base2}+=1;

		if (length($base1)>=$indelMaxLength){	$indelMaxLength = length($base1);	$indel=$base1;	};
		if (length($base2)>=$indelMaxLength){	$indelMaxLength = length($base2);	$indel=$base2;	};
	}else{
		$sum = $actPositionCounts{'A'} + $actPositionCounts{'T'} + $actPositionCounts{'C'} + $actPositionCounts{'G'} + $actPositionCounts{'N'} + $actPositionCounts{'X'};
		
		if ($sum == $desiredSum){	$type = "SNP";	}
		else{	$type = "INS";	}

		
		$n = $actPositionCounts{'N'};
		$dash = $actPositionCounts{'-'};	
#		$keys = join("\t", sort keys %actPositionCounts);		
		my $f = 1;
		
		foreach my $k (sort keys %actPositionCounts){
			if ($f == 1){
				$keys=$k;
				$values =$actPositionCounts{$k};
				$f=0;
			}else{
				$keys=$keys.",".$k;
				$values =$values.",".$actPositionCounts{$k};
			}
		}
		
						
		my $sig = "";
		if (exists($significant{$actPosition})){
			$sig = "SIG";	
		}else{
			$sig = "NON";
		}				
						
		print $actPosition."\t".$type."\t".$sig."\t".$indelMaxLength."\t".$keys."\t".$values."\n";

		$indelMaxLength=0;
		$indel="";		
		$type="";
		$actPosition=$pos;
		%actPositionCounts=();
		%actPositionCounts=('A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, 'N'=>0, 'X'=>0, '-'=>0);
		
		$actPositionCounts{$base1}+=1;
		$actPositionCounts{$base2}+=1;

		if (length($base1)>=$indelMaxLength){
			$indelMaxLength = length($base1);
			$indel=$base1;
		};
		if(length($base2)>=$indelMaxLength){
			$indelMaxLength = length($base2);
			$indel = $base2;
		};		

		
	}			
}
		

if ($sum == $desiredSum){	$type = "SNP";	}
else{	$type = "INS";	}

		
$n = $actPositionCounts{'N'};
$dash = $actPositionCounts{'-'};	
#$keys = join(",", keys %actPositionCounts);

my $f = 1;		
foreach my $k (sort keys %actPositionCounts){
	if ($f == 1){
		$keys=$k;
		$values =$actPositionCounts{$k};
		$f=0;
	}else{
		$keys=$keys.",".$k;
		$values =$values.",".$actPositionCounts{$k};
	}
}



my $sig = "";
if (exists($significant{$actPosition})){
	$sig = "SIG";	
}else{
	$sig = "NON";
}				

print $actPosition."\t".$type."\t".$sig."\t".$indelMaxLength."\t".$keys."\t".$values."\n";
#$indel."\t".$actPositionCounts{$indel}."\t".$n."\t".$dash."\t".$actPositionCounts{'A'}."\t".$actPositionCounts{'T'}."\t".$actPositionCounts{'C'}."\t".$actPositionCounts{'G'}."\t".$keys."\n";
		
close prettybaseFileHandle;


=pod

=head1 NAME 

summarize-prettybase.pl

=head1 SYNOPSIS

perl summarize-prettybase.pl --input fileName --list-of-significant listFileName

=head1 DESCRIPTION

The script reads an input prettybase file and prints out a summary of the file. The summary has form of a tab delimited file. For each position the following information is stored respectively in columns: 
	
	position, 
	type (INS for insertion, SNP for SNP), 
	significance (SIG for significant, NON for nonsignificant according to "--list-of-significant"), 
	max length of insertion on the position, 		
	a comma separated list of all alleles,
	a comma separated list of occurence counts for each allele (in the same order as the list of alleles).
	
=head1 OPTIONS 

=over 4

=item --input

An input file in prettybase format. It means, the file is tab delimited and consists of four columns: position, individual ID, allele1 and allele2. See, for example, http://www.pharmgat.org/Documentation/help/Prettybase. The option is mandatory.

=item --list-of-significant

A file that contains significant positions of an analysis we are interested in, optional. If the list is not provided, all the positions are considered to be nonsignificant.   

=item --number-of-individuals

Twice a number of diploid individuals used for the input prettybase file. The option is mandatory, default value is 192 (because of pigmentation, Bickel's data). 

=item --help, -h

Prints the help page.

=back

=head1 EXAMPLE

=over 4

=item INPUT file example.prettybase

	001469	id1	A	A
	001469	id2	T	T
	001469	id3	A	A
	001469	id4	A	A
	002637	id1	ACATCT	ACATCT
	002637	id2	ACATCT	ACATCT
	002637	id3	-	-
	002637	id4	N	N

=item INPUT file example.list

	2637

=item Used command

perl summarize-prettybase.pl --input example.prettybase --list-of-significant example.list --number-of-individuals 4 

=item OUTPUT

	1469	INS	NON	1	-,A,C,G,N,T,X	0,6,0,0,0,2,0
	2637	INS	SIG	6	-,A,ACATCT,C,G,N,T,X	2,0,4,0,0,2,0,0

=back 

=cut
