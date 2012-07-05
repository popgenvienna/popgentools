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


#  /Volumes/Volume_3/analysis/te/base/poly/Dmel-euchr-mc3-mq15-nr100-mincov10.txt  /Volumes/Volume_3/analysis/te/recombination/recombination-formated.txt 
my $poly_te_file    = shift;
my $rec_rate_file   = shift;
my $window          = shift; # 100_000
my $fixthreshold    = shift; # 0.95



my $pins=load_te_inserts($poly_te_file);
my $rr=load_recombination_rate($rec_rate_file);

my $teh=Utility::convert_tes_ins_hash($pins,$window);

while(my($chr,$temp)=each(%$teh))
{
    my $entries=@$temp;
    for my $i (0..$entries)
    {
        my $telist=$temp->[$i] || [];
        my $position=((2*$i+1)*$window)/2;
        my $recrate=$rr->{$chr}{$position} || 0;
        my $tecount=@$telist;
        my $fixed=0;
        foreach my $te (@$telist)
        {
            # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
            # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
            # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
            my $popfreq=$te->{popfreq};

            next unless $popfreq;
            $fixed++ if $popfreq > $fixthreshold;
        }
        my $fracfix=0;
        $fracfix=$fixed/$tecount if $tecount;
        print "$chr\t$position\t$tecount\t$fixed\t$fracfix\t$recrate\n";
    }
}

exit;

{
    use warnings;
    use strict;
    package Utility;
    
    sub convert_tes_ins_hash
    {
        my $telist=shift;
        my $winsize=shift;
        
        my $teinslist={};
        foreach my $te(@$telist)
        {
            # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
            # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
            # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
             my $chr=$te->{chr};
             my $pos=$te->{inspos};
             my $poskey=int($pos/$winsize);
             $teinslist->{$chr}[$poskey]=[] unless(exists($teinslist->{$chr}[$poskey]));
             push @{$teinslist->{$chr}[$poskey]},$te;             
        }
        return $teinslist;
    }
}