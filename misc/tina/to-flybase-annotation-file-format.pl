#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

my $fileName;
my $colPval;
my $HELP;

GetOptions(
	"input=s"=>\$fileName,
	"pvalue-column=i"=>\$colPval,
	'help'=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS"}) if $HELP;

if (!defined($fileName) or !defined($colPval)){die "Invalid arguments, use 'perl $0 --help'."}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

#print header of output file
my $header = "[general]\n".
			 "glyph = xyplot\n".
			 "graph_type=points\n".
			 "point_symbol=disc\n".
			 "point_radius=4\n".
			 "fgcolor = black\n".				 	
			 "bgcolor = black\n".
			 "height=150\n".
			 "min_score=0\n";
			 
			 
my @data;
my $max = 0;
my $printReference = 1;

			 
open fileHandle, "<", $fileName or die "Could not open file $fileName";
while (my $line = <fileHandle>){
	next if $line =~ m/^#/;
	
	chomp($line);

	my @parts = split "\t", $line;
	my $chr = $parts[0];
	my $pos = $parts[1];
	my $pval = $parts[$colPval-1];

	if ($printReference==1){
		print "reference=".$chr."\n";	
		$printReference=0;
	}

	my $logPval = -log10($pval);
	
	if ($logPval > $max){
		$max = $logPval;		
	}
	
	push @data, $chr."\t\"\"\t".$pos."..".$pos."\tscore=".$logPval."\n";			 

}
close fileHandle;
	
$header=$header."max_score=".$max."\n".
			 "label=0\n".	
			 "scale=left\n".
			 "key=-log_10(CMHtest p-value)\n";

print $header."\n";

foreach my $line (@data){
	print $line;	
}


=pod

=head1 NAME

to-flybase-annotation-file-format.pl

=head1 SYNOPSIS

perl to-flybase-annotation-file-format.pl --input CMHtestFile --pvalue-column 9 

=head1 DESCRIPTION 

The script takes an input file with CMH test p-values and convert it to a file for flybase. After import of the file to flybase,
manhattan plot of a small region can be displayed together with other flybase features.

=head1 OPTIONS

=over 4 

=item --input 

A CMHtest output file.

=item --pvalue-column

A number of column that contains CMH test p-value. Columns in the file are numbered form 1. 

=item --help, -h

Prints the help page.

=back
