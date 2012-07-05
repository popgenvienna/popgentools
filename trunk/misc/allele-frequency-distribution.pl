use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use Pileup;


my $input;
my $output;
my $mincount=2;
my $minqual=20;
my $mincov=4;
my $maxcov=1000000;
my $bincount=20;

my $help=0;
GetOptions(
           "input=s"        => \$input,
           "output=s"       => \$output,
           "min-count=i"    => \$mincount,
           "min-qual=i"     => \$minqual,
           "min-cov=i"      => \$mincov,
           "max-cov=i"      => \$maxcov,
           "bin-number=i"   => \$bincount,
           "help"           => \$help
        ) or pod2usage(-msg=>"Invalid arguments",-verbose=>1);


pod2usage(-verbose=>2) if $help;
pod2usage(-msg=>"Could not find input file file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Output file must be specified",-verbose=>1) unless $output;

my $bins=_getbins($bincount);

open my $ifh,"<",$input or die "Could not open input file";
open my $ofh,">",$output or die "Could not open output file";

# qualencoding,mincount,mincov,maxcov,minqual
my $pp=get_pileup_parser("illumina",$mincount,$mincov,$maxcov,$minqual);

my $countNoRefAllele=0;
while(my $l=<$ifh>)
{
    chomp $l;
    # line, mincount, mincoverage, minqual, maxcoverage
    my $p=$pp->($l);
    next unless $p->{ispuresnp};
    my $rc=$p->{refc};
    next if $rc eq "N";
    
    if($p->{$rc}==0)
    {
        $countNoRefAllele++;
        next;
    }

    my $refcfreq=_getrefcharfreq($p);
    
    foreach my $b (@$bins)
    {
        if($refcfreq>=$b->{lower} and $refcfreq< $b->{upper})
        {
            $b->{count}++;
        }
    }
    
}
    print $ofh "SNPs with zero reference allele count: $countNoRefAllele\n";
    # print the different bins
    foreach my $b (@$bins)
    {
        print $ofh "$b->{lower}\t$b->{upper}\t$b->{count}\n";
    }

exit;

sub _getrefcharfreq
{
    my $snp=shift;
    my $rc=$snp->{refc};
    
    my $refcount=$snp->{$rc};
    my $reffreq=$refcount/$snp->{eucov};
    return $reffreq;
}

sub _getbins
{
    my $bin=shift;
    my $step=1/$bin;
    
    my $bins=[];
    my $lowerbound=0;
    my $upperbound=$step;
    foreach(1..$bin)
    {
        push @$bins,{
            lower=>$lowerbound,
            upper=>$upperbound,
            count=>0
            };
        $lowerbound+=$step;
        $upperbound+=$step;
    }
    #$bins->[-1]{upper}+=0.000001;
    return $bins;
}

=head1 NAME

perl allele-frequency-distribution.pl - A script which calculates the allele frequency spectrum

=head1 SYNOPSIS

perl allele-frequency-distribution.pl --input input.pileup --output distribution.txt

=head1 OPTIONS

=over 4

=item B<--input>

The input file;  a pileup-file;  Mandatory.

=item B<--output>

The allele frequency distribution; Mandatory

=item B<--min-count>

The minimum count of the minor allele; default=2

=item B<--min-cov>

The minimum coverage; default=4

=item B<--max-cov>

The maximum coverage; default=1000000

=item B<--min-qual>

The minimum quality; default=20

=item B<--bin-number>

The number of bins to use for this analysis; default=20

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

a psl file

=head2 Output

a gtf file
  
=cut