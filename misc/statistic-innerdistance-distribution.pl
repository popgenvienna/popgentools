#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use SAMPairReader;


my $input;
my $minMapQual=0;
my $help=0;

my $maxdist=1000;

    
GetOptions(
           "input=s"        => \$input,
           "max-dist=i"     => \$maxdist,
           "help"           => \$help   );

pod2usage(-verbose=>2) if $help;

my $sr=SAMPairReader->new($input,10000000000);

my $idh={};
while(my $pair=$sr->next())
{
    my($r1,$r2) = @$pair;
    # readid, flag, chr, start, mq, cigar, chrmate, posmate, distance, seq, qual, appendix,end, start_s, end_s
    next unless $r1->{chr} eq $r2->{chr};
    next if(($r1->{flag} & 0x0010) == ($r2->{flag} & 0x0010));
    
    my $innerdist=$r2->{start_s}-$r1->{end_s};
    
    next if(abs($innerdist)> $maxdist);
    $idh->{$innerdist}++;
}


my $toprintar=[];
while(my($dist,$count)=each(%$idh)){
    push @$toprintar,{
        id=>$dist,
        count=>$count
    };
}

$toprintar=[sort {$a->{id}<=>$b->{id}} @$toprintar];

foreach my $t (@$toprintar)
{
    print "$t->{id}\t$t->{count}\n";
}

exit;

=head1 NAME

perl statistic-innerdistance-distribution.pl - A script calculating the distribution of the inner distance given a sam file

=head1 SYNOPSIS

perl statistic-innerdistance-distribution.pl --max-dist 1000 --input pe_sam.sam

=cut