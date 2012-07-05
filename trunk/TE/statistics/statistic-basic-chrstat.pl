use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;

my $file=shift;
my $minfix=shift;

unless($minfix)
{
    print "$0 filename minfix\n";
    exit;
}



my $teins=load_te_inserts($file);

Utility::statistics($teins,$minfix);


exit;




{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    sub statistics
    {
        # print for each chromosome the number of transposable elements;

        # 
        my $teins=shift;
        my $minfix=shift;
        
        my $tedistri={};

        # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
        foreach my $te (@$teins)
        {
            my $pos=$te->{inspos};
            my $popfreq=$te->{popfreq};
            my $fix=$popfreq> $minfix?1:0;
            # Position key
            my $chr=$te->{chr};
            $tedistri->{$chr}{all}++;
            $tedistri->{$chr}{fix}++ if $fix;
        }
        
        
        while(my($chr,$temp)=each(%$tedistri))
        {
            my $chrdist_all=$temp->{all};
            my $chrdist_fix=$temp->{fix};
            print "$chr\t$chrdist_all\t$chrdist_fix\n";
        }
        
    }
}





