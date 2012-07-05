{
    package Recombination;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(load_recombination_rate);
    
    
    sub load_recombination_rate
    {
        my $file=shift;
        open my $ifh, "<", $file or die "could not open file";
        
        my $rr={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            my($chr,$pos,$recr)=split /\t/,$l;
            $rr->{$chr}{$pos}=$recr;
        }
        return $rr;
    }
    
}