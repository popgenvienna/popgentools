#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;


my $popgentoolpath="/Users/robertkofler/dev/PopGenTools/";
my $popoolationpath="/Users/robertkofler/dev/popoolation/";

my @dirs=qw{basic-pipeline Modules Modules/Test};
my @files=qw{calculate-dxy.pl VarSliding2Wiggle.pl VarSliding2Flybase.pl mauve-parser.pl Variance-at-position.pl Variance-sliding.pl Visualise-output.pl
 basic-pipeline/filter-pileup-by-gtf.pl basic-pipeline/identify-genomic-indel-regions.pl basic-pipeline/mask-sam-indelregions.pl basic-pipeline/trim-fastq.pl basic-pipeline/convert-fastq.pl basic-pipeline/subsample-pileup.pl
 Modules/BasicUtility.pm Modules/Pileup.pm Modules/Test.pm Modules/VarianceExactCorrection.pm Modules/VarianceUncorrected.pm Modules/VarMath.pm Modules/FormatSNP.pm
 Modules/Test/PileupParser.pm Modules/Test/Variance.pm};


foreach my $d (@dirs)
{
    my $fuld=$popoolationpath.$d;
    unless(-e $fuld)
    {
        mkdir($fuld) or die "could not create directory";
        print "Created directory\n";
    }
}

foreach my $f (@files)
{
    my $from=$popgentoolpath.$f;
    my $to=$popoolationpath.$f;
    die "$from does not exist" unless -e $from;
    copy($from, $to) or die "File $from cannot be copied to $to";
    print "Copied file from $from to $to\n";
    
}
