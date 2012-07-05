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
use Getopt::Long;


my $teinsertfile=shift;
my $lengthhashfile=shift;


my $lengannotator=Utility::get_telength_annotator($lengthhashfile);




my $polyins=load_te_inserts($teinsertfile);
foreach my $ins (@$polyins)
{
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap

    my $teid = $ins->{teid};
    my $chr  =  $ins->{chr};
    my $pos  = int($ins->{inspos});
    
    my $fbid = $ins->{fbid};
    next unless $fbid;
    next unless $lengannotator->($fbid);

    my $leng=$lengannotator->($fbid);
    print  $ins->tostring()."\t$leng\n";
}

exit;




{
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use List::Util qw[min max];
    use Recombination;
    package Utility;
    
    

    
    sub get_telength_annotator
    {
        my $telengthfile=shift;
        my $famlengh={};
        open my $ifh, "<",$telengthfile or die "Could not open input file";
        while(my $l=<$ifh>)
        {
            chomp $l;
            my($fam,$leng)=split /\t/,$l;
            $famlengh->{$fam}=$leng;
        }
        
        return sub{
            my $queryfam=shift;
            die "no entry for family $queryfam" unless(exists($famlengh->{$queryfam}));
            return $famlengh->{$queryfam};
        }
    }
    
    
    

}