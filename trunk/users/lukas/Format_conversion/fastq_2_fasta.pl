#!/usr/bin/perl
#
# This program converts Sanger formatted fastq files to fasta files
# The Sanger format encodes quality scores 0-93 as ASCII values 33-126
# It is not compatable with the Illumina fastq format.
#
# Revised on 11/15/09 to fix bug in how masking coordinates
# 

use Getopt::Std;

getopt("qoim");


$usage = <<END;
	-i input fastq file
	-o output fasta file
	-q [0] quality score cutoff 
	-m mask file (line chr start stop)

Note: For optional mask file the chromosome begins at  
      position 1 and the range [start stop] is masked.
END

die $usage unless (-e $opt_i && $opt_o);
open IN, "$opt_i" || die "Could not open $opt_i\n";
open OUT, ">$opt_o" || die "Could not open $opt_o\n";

if ($opt_q > 0 && $opt_q <= 93) {
	$opt_q += 33;
}
else {
	$opt_q = 0;
}

if ( -e $opt_m ) {
	open MASK, $opt_m;
	while (<MASK>) {
		next if /^#/;
		if ( /(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/ ) {
			push @{ $masks{"$1\t$2"} }, [$3, $4];
		}
	}
}

while(<IN>) {
#        if ( /^@(\S+)\s+(\S+)/ )
	if ( /^@(\w+)/ ) {
#		$id = "$1\t$2";
	        $id = "$1";
		$seq = <IN>;
		chomp $seq;
		<IN>;
		$qual = <IN>;
		chomp $qual;
		#assert( length($seq) == length($qual) );
		if ( $opt_q ) {
			for ($i = 0; $i < length($qual); $i++) {
				if ( ord(substr($qual, $i, 1)) < $opt_q ) {
					substr($seq, $i, 1) = "N";
				}
			} 
		}
		if ( $opt_m ) {
			foreach ( @{ $masks{$id} } ) {
				substr($seq, $_->[0]-1, $_->[1]-$_->[0]+1) = "N" x ($_->[1]-$_->[0]+1);
			}
		}
		print OUT ">$id\n$seq\n";
	}
}
