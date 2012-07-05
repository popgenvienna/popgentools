#/usr/bin/perl 

use warnings;
use strict;

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

my $HELP;

GetOptions(
  'help'=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|EXAMPLE"}) if $HELP;

my %hash;

while (my $line = <STDIN>){
	chomp($line);
	my ($chr, $junk, $junk2, $start, $end) = split "\t", $line;

	if(! defined($hash{$chr})){
		$hash{$chr}{chr}=$chr;
		$hash{$chr}{min}=$start;
	#	$hash{$chr}{min}=1;
		$hash{$chr}{max}=$end;
	}else{
		if ($hash{$chr}{min}>$start){
			$hash{$chr}{min}=$start
		}
		if ($hash{$chr}{max}<$end){
			$hash{$chr}{max}=$end	
		}	
	}
}

foreach my $chromosome (keys %hash){
	my $string = $chromosome."\t"."create-chromosome-arms.pl"."\t"."chromosome_arm"."\t".$hash{$chromosome}{min}."\t".$hash{$chromosome}{max}."\t.\t.\t.\n";
	print $string;
}

=pod

=head1 NAME 

create-chromosome-arms.pl

=head1 SYNOPSIS

perl create-chromosome-arms.pl < inputGFF/GTFfile 

=head1 DESCRIPTION

The script takes an input file in GFF/GTF format and creates chromosome_arm ferature line for each of chromosome arms in the input file.

=head1 EXAMPLE

=over 4

=item INPUT gff file example.gff

  2L	flybase	exon	15	36	.	.	.
  2L	flybase	intron	38	50	.	.	.
  2L	flybase	exon	55	70	.	.	.
  
=item used command

perl create-chromosome-arms.pl < example.gff

=item OUTPUT
  
  2L	create-chromosome-arms.pl	chromosome_arm	15	70	.	.	.

=back

=cut  
