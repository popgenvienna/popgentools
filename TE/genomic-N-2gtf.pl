#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use Getopt::Long;
use Pod::Usage;
use Utility;

my $input;
my $help=0;

GetOptions(
    "input=s"       =>\$input,
    "help"          =>\$help
);
pod2usage(-verbose=>2) if $help;

my $fastr = get_fasta_reader($input);
my $counter=1;
while(my $fast=$fastr->())
{
    my $header=$fast->{head};
    $header=~s/\s+.*^//;
    my $seq=$fast->{seq};
    while($seq=~m/(N+)/gi)
    {
        # ANNNNAATCG
        my $end=pos($seq);
        my $len=length($1);
        my $start=$end-$len+1;
        #2R	maker	CDS	19450390	19450625	.	-	0	ID
        print "$header\tpolyNsearch\tpolyN\t$start\t$end\t.\t+\t.\tgene_id \"poly_N_$counter\";\n";
        $counter++;
    }
}


=head1 NAME

perl genomic-N-2gtf.pl - Identifies the position of poly-N stretches and outputs these positions in a gtf file

=head1 SYNOPSIS

 perl genome-N-2gtf.pl --input whole_genome.fa 

=head1 OPTIONS

=over 4

=item B<--input>

the reference sequence (after repeat-masking) in fasta format. Mandatory

=item B<--output>

The output file. Mandatory

=cut