use strict;
use warnings;
use Getopt::Long;
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use List::Util qw[min max];
use TEInsertUtility;


my $contigsonly=[qw/2L 2R 3L 3R 4 X/];
my $contighash= {map {($_,1) } @$contigsonly};


my $teinsert_file;
my $variance_file;
my $stepwinsize;
my $range=3;
my $varthreshold;

GetOptions(
    "insert-file=s"	    =>\$teinsert_file,
    "variance-file=s"       =>\$variance_file,
    "step-window-size=i"    =>\$stepwinsize,
    "threshold=f"           =>\$varthreshold
) or die "Wrong parameters";


my $teinserts=load_te_inserts($teinsert_file);
my $variance=Utility::load_variance($variance_file, $stepwinsize);


foreach my $tei (@$teinserts)
{
    
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
    my $chr=$tei->{chr};
    my $pos=$tei->{inspos};
    
    my $relpos=int($pos/$stepwinsize);
    my @positions=(($relpos-$range)..($relpos+3));
    while($positions[0]<1)
    {
        shift @positions;
    }
    my $varar=[];
    
    my $toprint=0;
    foreach my $i (@positions)
    {
        my $var=$variance->{$chr}[$i];
        $var="na" unless(defined($var));
        if(($i==$relpos or $i-1==$relpos or $i+1==$relpos) and $var ne "na")
        {
            $toprint=1 if $var <=$varthreshold;
        }
        $var= "'$var'" if $i == $relpos;
        push @$varar,$var;
    }
    next unless $toprint;
    my $varstring =join (" ",@$varar);
    
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
    my $fbid=$tei->{fbid}||"-";
    print "$tei->{chr}\t$tei->{inspos}\t$tei->{sitesupport}\t$tei->{teid}\t$tei->{popfreq}\t$tei->{order}\t$fbid\t$varstring\n";
    
    
}



my $bla=0;

exit;




{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    sub load_variance
    {
        my $varfile     =shift;
        my $stepwinsize =shift;
        
        open my $ifh,"<",$varfile or die "Could not open variance file"; 
        my $varhash={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            # 2L	633150	1	1.000	0.003489767
            # 2L	633250	1	1.000	0.005335821
            # 2L	633350	2	1.000	0.006082731
            my($chr,$pos,$snps,$covfrac,$var)=split /\t/,$l;
            my $relpos=int($pos/$stepwinsize);
            die "" if(exists($varhash->{$chr}[$relpos]));
            $varhash->{$chr}[$relpos]=$var;
        }
        return $varhash;
    }
    

    

    
 
    

    
 
}





