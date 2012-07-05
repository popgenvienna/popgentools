use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;

my $file=shift;
my $family=shift;
my $windowsize=shift;
my $minfix=shift;

unless($minfix)
{
    print "$0 filename windowsize minfix\n";
    exit;
}



my $teins=load_te_inserts($file);

Utility::print_polynonpoly_correlation($teins,$windowsize,$minfix);


exit;




{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    sub print_polynonpoly_correlation
    {
        # print for each chromosome the number of transposable elements;

        # 
        my $teins=shift;
        my $windowsize=shift;
        my $minfix=shift;
        
        my $tedistri={};

        # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
        foreach my $te (@$teins)
        {
            my $tfam=$te->{teid};
            next unless $tfam eq $family;
            my $pos=$te->{inspos};
            my $popfreq=$te->{popfreq};
            my $fix=$popfreq> $minfix?1:0;
            # Position key
            my $poskey=int($pos/$windowsize);
            my $chr=$te->{chr};
            $tedistri->{$chr}{all}[$poskey]++;
            $tedistri->{$chr}{fix}[$poskey]++ if $fix;
        }
        
        

        
        while(my($chr,$temp)=each(%$tedistri))
        {
            my $chrdist_all=$temp->{all};
            my $chrdist_fix=$temp->{fix};
            my $leng=max(scalar(@$chrdist_all),scalar(@$chrdist_fix));
            foreach my $i (0..($leng-1))
            {
                my $pos=(2*$i+1)*$windowsize/2;
                
                my $cdall   = $chrdist_all->[$i];
                my $cdfix   = $chrdist_fix->[$i];
                
                $cdall ||=0;
                $cdfix ||=0;
                
                print "$chr\t$pos\t$cdall\t$cdfix\n";
            }
        }
        
    }
}





