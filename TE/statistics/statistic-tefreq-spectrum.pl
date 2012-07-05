use strict;
use warnings;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;

my $classes=10;


my $file=shift;

die "$0 filename" unless -e $file;

my $teh={};
open my $ifh, "<", $file or die "Could not open input file";

my $teins=load_te_inserts($file);

foreach my $te (@$teins)
{
            # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
        # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
        # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
    my $teid=$te->{teid};
    my $freq=$te->{popfreq};
    my $cat=sprintf("%.1f",int($freq*$classes)/$classes);
    $teh->{$teid}{$cat}++;
}

while(my($te,$temp)=each(%$teh))
{
    my $abs=0;
    my $counter=0;
    my @freqs=();
    for my $i (0..10)
    {
        my $cat=sprintf("%.1f",$i/$classes);
        my $val=$temp->{$cat};
        $val||=0;
        $counter+=$val;
        $abs+=($val*$cat);
        push @freqs,$val;
    }
    my $spectr=$abs/$counter;
    my $freqs=join(",",@freqs);
    $spectr=sprintf("%.3f",$spectr);
    print "$te\t$counter\t$spectr\t$freqs\n";
}