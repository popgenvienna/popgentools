#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use ParseSam;
use PESamReader;


my $input=shift;
my $tosample=shift;
my $insample=shift;


die "usage: $0 input tehierfile tosample insample" unless $input;

my $sr=PESamReader->new($input);

my $ratio;

while(my $s=$sr->next())
{
        $ratio=$tosample/$insample;
        my $rand=rand();
        if($rand<$ratio)
        {
            print $s->[0]."\n";
            print $s->[1]."\n";

            $tosample--;
        }
    
    $insample--; #pe fragments every fragments only counts for half
}

$sr->close()









