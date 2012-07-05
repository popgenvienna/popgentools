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

my $famhash={};

# distribute each family member into an array
foreach my $te (@$teins)
{
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
    my $teid=$te->{teid};
    $famhash->{$teid}||=[];
    push @{$famhash->{$teid}},$te
  
}

print "\tfamily\tcount\tavfreq\tavpi\tsumfreq\tsumpi\n";
while(my($fam,$elements)=each(%$famhash))
{
    my $count=@$elements;
    my $freqsum=0;
    my $pisum=0;
    
    
    for my $e (@$elements)
    {
        my $freq=$e->{popfreq};
        $freqsum+=$freq;
        my $pi=1-$freq**2-(1-$freq)**2;
        $pisum+=$pi;
    }
    
    my $averagefreq=sprintf("%.3f",$freqsum/$count);
    my $averagepi=sprintf("%.3f",$pisum/$count);
    print "$fam\t$count\t$averagefreq\t$averagepi\t$freqsum\t$pisum\n";
}