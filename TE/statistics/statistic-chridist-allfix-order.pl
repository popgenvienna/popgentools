use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;


my $file=shift;
my $windowsize=shift;
my $minfix=shift;
print "$0 filename windowsize minfix\n" unless $minfix;

my $teins=load_te_inserts($file);


Utility::print_te_polynonpoly($teins,$windowsize,$minfix);


exit;




{
    package Utility;
    use strict;
    use warnings;
    
    
    sub print_te_polynonpoly
    {

    
        my $teins=shift;
        my $windowsize=shift;
        my $minfix=shift;
        
        my $tedistri={};
        my $maxkey=0;

        foreach my $te (@$teins)
        {
            # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
            
            # choose any identifier here (order or teid, family)
            my $teid=$te->{order};
            
            my $position=$te->{inspos};
            my $popfreq=$te->{popfreq};            
            my $fix=$popfreq>$minfix?1:0;

            # Position key
            my $poskey=int($position/$windowsize);
            $maxkey =$poskey if $poskey > $maxkey;
            
            my $chr=$te->{chr};
            $tedistri->{$teid}{$chr}{all}[$poskey]++;
            $tedistri->{$teid}{$chr}{fix}[$poskey]++ if $fix;
            
        }
        
        # print header
        print "teid\tchr\tcat";
        for my $i (1..$maxkey)
        {
            my $start=$windowsize*($i-1);
            my $end=$windowsize*$i;
            my $pos=($end+$start)/2;
            print "\t$pos";
        }
        print "\n";
        
        
        while(my($teid,$temp2)=each(%$tedistri))
        {
            print "$teid\n";
            while(my($chr,$temp)=each(%$temp2))
            {
                foreach my $cat (qw/all fix/)
                {
                    my $chrdist=$temp->{$cat};
                    next unless $chrdist;
                    my $leng=@$chrdist;
                    my @vals=();
                    foreach my $i (0..($leng-1))
                    {
                        my $c=$chrdist->[$i];
                        $c||=0;
                        push @vals,$c;
                    }
                    @vals=map{sprintf("%01d",$_) } @vals;
                    my $chrstring=join("\t",@vals);
                    print "$teid\t$chr\t$cat\t$chrstring\n";
                }
                print "\n";
            }
            print "\n";
        }
    }
}





