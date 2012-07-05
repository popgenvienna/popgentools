use strict;
use warnings;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use List::Util qw[min max];
use TEInsertUtility;
use Recombination;

my $classes=10;
my $windowsize=100000;

my $outputprefix    = shift;
my $polytes         = shift;
my $rec_rate_file   = shift;
my $recthreshold    = shift;


my $output_low=$outputprefix."_low";
my $output_high=$outputprefix."_high";

open my $ofhl,">",$output_low or die "Could not open output file";
open my $ofhh,">",$output_high or die "Could not open output file";

my $teh={};

my $polyins=load_te_inserts($polytes);
my $rr=load_recombination_rate($rec_rate_file);


foreach my $ins (@$polyins)
{
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap

    my $teid = $ins->{teid};
    my $pos  = $ins->{inspos};
    my $index =int($pos/$windowsize);
    my $reckey = ((2*$index+1)*$windowsize)/2;
    my $recatpos=$rr->{$ins->{chr}}{$reckey} || 0;
    
    if($recatpos>$recthreshold)
    {
        print $ofhh $ins->tostring()."\n"; 
    }
    else
    {
        print $ofhl $ins->tostring()."\n";
    }

}

exit;




{
    use strict;
    use warnings;
    package Utility;
    
    sub get_median
    {
        my $rr=shift;
        my @rar=();
        while(my($chr,$temp)=each(%$rr))
        {
            while(my($pos,$recr)=each(%$temp))
            {
                push @rar,$recr;
            }
        }
        @rar = sort {$a<=>$b} @rar;
        my $index = int(scalar(@rar/2));
        
        return $rar[$index];
    }
    
    sub write_output
    {
        my $ins=shift;
        #2L      F       INE-1   20922966        20923040        0       40      0
        #2L      R       INE-1   20923660        20923734        0       19      0
        # chr, insdir, teid, start, end, count_pre, count_abs
        my $cpre=$ins->{count_pre};
        my $cabs=$ins->{count_abs};
        my $pi=$ins->{pi} ||0;
        return "$ins->{chr}\t$ins->{insdir}\t$ins->{teid}\t$ins->{start}\t$ins->{end}\t$pi\t$cpre\t$cabs\n";   
    }
}