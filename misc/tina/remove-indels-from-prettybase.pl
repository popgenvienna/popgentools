#!/usr/bin/perl 

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $prettybase;
my $desiredSum = 192;
my $HELP;

GetOptions(
	"input=s"=>\$prettybase,
	"number-of-individuals=i"=>\$desiredSum,
	"help"=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;
 
my @onePositionData=();
my $actPosition;
my %actPositionCounts=('A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, 'N'=>0, 'X'=>0, '-'=>0); 

my $sum;
my $noOfPositions = 0;
my $noOfIndels = 0;

my $indelMaxLength=0;
my $indel;

 
open prettybaseFileHandle, "<", $prettybase;

my $line = <prettybaseFileHandle>;
chomp($line);
my ($pos, $flyID, $base1, $base2)	 = split "\t", $line;
$pos =~ s/^0*//;

$actPosition = $pos;
push (@onePositionData, $pos."\t".$flyID."\t".$base1."\t".$base2);
$actPositionCounts{$base1}+=1;
$actPositionCounts{$base2}+=1;


my $n=0;
my $dash=0;
my $keys;
my $values;


print STDERR "#indelPosition\tmaxLengthOfTheInsertionOrDeletion\talleles\tcounts\n";


while ($line = <prettybaseFileHandle>){
	chomp($line);
	my ($pos, $flyID, $base1, $base2)	 = split "\t", $line;
	$pos =~ s/^0*//;

	if ($pos == $actPosition){
		$actPosition = $pos;
		push (@onePositionData, $pos."\t".$flyID."\t".$base1."\t".$base2);
		$actPositionCounts{$base1}+=1;
		$actPositionCounts{$base2}+=1;
		
		if (length($base1)>=$indelMaxLength){	$indelMaxLength = length($base1);	};
		if(length($base2)>=$indelMaxLength){	$indelMaxLength = length($base2);	};
		
	}else{
		
		$sum = $actPositionCounts{'A'} + $actPositionCounts{'T'} + $actPositionCounts{'C'} + $actPositionCounts{'G'} + $actPositionCounts{'N'} + $actPositionCounts{'X'};
				
		if ($sum == $desiredSum){
			print join("\n", @onePositionData);	
			print "\n";
			$noOfPositions+=1;
		}else{

			$n = $actPositionCounts{'N'};
			$dash = $actPositionCounts{'-'};	
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

						
			print STDERR $actPosition."\t".$indelMaxLength."\t".$keys."\t".$values."\n";
			$noOfIndels+=1;
		}

		$indelMaxLength=0;
		$indel="";		
		$actPosition=$pos;
		%actPositionCounts=();
		%actPositionCounts=('A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, 'N'=>0, 'X'=>0, '-'=>0);
		@onePositionData=();
		
		push (@onePositionData, $pos."\t".$flyID."\t".$base1."\t".$base2);
		$actPositionCounts{$base1}+=1;
		$actPositionCounts{$base2}+=1;

		if (length($base1)>=$indelMaxLength){	$indelMaxLength = length($base1);	};
		if(length($base2)>=$indelMaxLength){	$indelMaxLength = length($base2);	};		
	}
}

$sum = $actPositionCounts{'A'} + $actPositionCounts{'T'} + $actPositionCounts{'C'} + $actPositionCounts{'G'} + $actPositionCounts{'N'} + $actPositionCounts{'X'};

if ($sum == $desiredSum){
	print join("\n", @onePositionData);	
	print "\n";
	$noOfPositions+=1;
}else{
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
	print STDERR $actPosition."\t".$indelMaxLength."\t".$keys."\t".$values."\n";
	$noOfIndels+=1;
}

print "#".$noOfPositions, "\n";
print STDERR "#".$noOfIndels,"\n";

close prettybaseFileHandle;


=pod

=head1 NAME

remove-indels-from-prettybase.pl

=head1 SYNOPSIS

perl remove-indels-from-prettybase.pl --input fileName --number-of-individuals number > outFileName 2> errFileName

=head1 DESCRIPTION

The script takes input prettybase file and prints out prettybase records only for the SNP positions (no insertions, no deletions). This lines are printed to standard output, the script also prints some additional information to standard error output.
Summary information about each nonSNP position is printed to error output. At the end of both outputs a number of positions is printed in a form of commented line (starting with '#'). 

=head1 OPTIONS

=over 4

=item --input

An input file in prettybase format. It means, the file is tab delimited and consists of four columns: position, individual ID, allele1 and allele2. See, for example, http://www.pharmgat.org/Documentation/help/Prettybase for details. The option is mandatory.

=item --number-of-individuals

Twice a number of diploid individuals used for the input prettybase file. The option is mandatory, default value is 192 (historically, because of Bickel's pigmentation data). 

=item --help, -h 

=back

Prints the help page.

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

=item Used command

perl remove-indels-from-prettybase.pl --input example.prettybase --number-of-individuals 8 > output.prettybase 2> err.indels

=item OUTPUT file output.prettybase

	1469	id1	A	A
	1469	id2	T	T
	1469	id3	A	A
	1469	id4	A	A
	#1

=item OUTPUT file err.indels

	#indelPosition	maxLengthOfTheInsertionOrDeletion	alleles	counts
	2637	6	-,A,ACATCT,C,G,N,T,X	2,0,4,0,0,2,0,0
	#1

=back

=cut
