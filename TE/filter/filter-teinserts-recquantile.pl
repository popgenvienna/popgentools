use strict;
use warnings;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;
use Recombination;

#my $classes=10
my $windowsize=100000;

my $outputprefix    = shift;
my $teinserts       = shift;
my $rec_rate_file   = shift;
my $recquantiles    = shift;

die "Recombination quantiles are not set $recquantiles" unless $recquantiles;

my $rq      = Utility::get_recombination_quantiles($outputprefix,$recquantiles);
my $polyins = load_te_inserts($teinserts);
my $rr      = load_recombination_rate($rec_rate_file);

foreach my $tei (@$polyins)
{
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
    my $teid = $tei->{teid};
    my $pos  = $tei->{inspos};
    my $index =int($pos/$windowsize);
    my $reckey = ((2*$index+1)*$windowsize)/2;
    my $recatpos=$rr->{$tei->{chr}}{$reckey} || 0;
    
    foreach my $rec (@$rq)
    {
        if($recatpos>=$rec->{lower} and $recatpos < $rec->{upper})
        {
            my $fh=$rec->{fh};
            print $fh $tei->tostring()."\n";
        }
    }
    

}

exit;




{
    use strict;
    use warnings;
    package Utility;
    
    sub get_recombination_quantiles
    {
        my $outputprefix = shift;
        my $recquantiles = shift;
        my @rar=split /\s/,$recquantiles;
        
        my $toret=[];
        
        for my $i(0..($#rar-1))
        {
            my $lower=$rar[$i];
            my $upper=$rar[$i+1];
            
            my $filename=$outputprefix."_$lower"."_$upper";
            
            open my $ofh,">",$filename or die "could not open output file $filename";
            push @$toret,{
                filename=>$filename,
                fh=>$ofh,
                lower=>$lower,
                upper=>$upper,
            };
        }
        return $toret;
    }
    

    
    sub write_output
    {
        my $ins=shift;
        #2L      F       INE-1   20922966        20923040        0       40      0
        #2L      R       INE-1   20923660        20923734        0       19      0
        # chr, insdir, teid, start, end, count_pre, count_abs
        my $cpre=$ins->{count_pre} ;
        my $cabs=$ins->{count_abs} ;
        my $pi=$ins->{pi} ||0;
        return "$ins->{chr}\t$ins->{insdir}\t$ins->{teid}\t$ins->{start}\t$ins->{end}\t$pi\t$cpre\t$cabs\n";   
    }
}