#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;


my $chrinfo={
    X       =>{leng=>22422827, cent_pos=>"end"},
    '2L'    =>{leng=>23011544, cent_pos=>"end"},
    '2R'    =>{leng=>21146708, cent_pos=>"start"},
    '3L'    =>{leng=>24543557, cent_pos=>"end"},
    '3R'    =>{leng=>27905053, cent_pos=>"start"}
    };

my $teinsertfile=shift;
my $recombinationrate=shift;
my $recthreshold=shift;
my $outputprefix=shift;

my $fhhash={};
open my $tfh,">",$outputprefix."_telo";
$fhhash->{telo}=$tfh;
open my $cfh,">",$outputprefix."_cent";
$fhhash->{cent}=$cfh;
open my $efh,">",$outputprefix."_euchr";
$fhhash->{euchr}=$efh;

my $reghash=Utility::define_regions($chrinfo,$recombinationrate,$recthreshold);
my $teins=load_te_inserts($teinsertfile);

foreach my $tei (@$teins)
{
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
    
    my $chr=$tei->{chr};
    next unless(exists($chrinfo->{$chr}));
    my $pos=$tei->{inspos};
    
    #get region
    my $region;
    foreach my $candreg (@{$reghash->{$chr}})
    {
        if($pos>$candreg->{start} and $pos<=$candreg->{end})
        {
            $region=$candreg;
        }
    }
    die ""unless $region;
    

    my $feat=$region->{feat};
    my $fh=$fhhash->{$feat};
    print $fh $tei->tostring()."\n";
}





{
    package Utility;
    use strict;
    use warnings;
    
    sub load_recombination_rate
    {
        my $recrate=shift;
        open my $ifh, "<",$recrate or die "Could not open recombination rate";
        
        my $rh={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            my($chr,$pos,$rr)=split /\t/,$l;
            $rh->{$chr}||=[];
            push @{$rh->{$chr}},{
                pos=>$pos,
                rr=>$rr
            };
        }
        #2L	50000	0.00
        #2L	150000	0.00
        #2L	250000	0.00
        return $rh;
    }
    
    

    sub define_regions
    {
        my $chrinfo=shift;
        my $recombinationfile=shift;
        my $recthreshold=shift;
        
        my $rh=load_recombination_rate($recombinationfile);
        
        my $reghash={};
        while(my($chr,$e)=each(%$chrinfo))
        {
            my $cent_pos=$e->{cent_pos};
            my $chrrec=$rh->{$chr};
            my $leng=$e->{leng};
            
            my $leftpos=0;
            my $rightpos=$e->{leng};
            
            while(@$chrrec)
            {
                my $ca=shift @$chrrec;
                my $trec=$ca->{rr};
                if($trec < $recthreshold)
                {
                  $leftpos=$ca->{pos};  
                }
                else
                {
                    unshift @$chrrec,$ca; 
                    last;
                }
            }
            while(@$chrrec)
            {
                my $ca=pop @$chrrec;
                my $trec=$ca->{rr};
                if($trec < $recthreshold)
                {
                    $rightpos=$ca->{pos};
                }
                else
                {
                    push @$chrrec,$ca;   
                    last;
                }
            }
            
            my $feat=[];
            push @$feat,
            {
                feat=>"euchr",
                start=>$leftpos,
                end=>$rightpos
            };
            if($cent_pos eq "end")
            {
                push @$feat,
                {
                    feat=>"cent",
                    start=>$rightpos,
                    end=>$leng
                };
                push @$feat,
                {
                    feat=>"telo",
                    start=>0,
                    end=>$leftpos
                };
                
            }
            elsif($cent_pos eq "start")
            {
                push @$feat,
                {
                    feat=>"telo",
                    start=>$rightpos,
                    end=>$leng
                };
                push @$feat,
                {
                    feat=>"cent",
                    start=>0,
                    end=>$leftpos
                };
            }
            $_->{count_all}=0 foreach @$feat;
            $_->{count_fix}=0 foreach @$feat;
            $_->{leng}=$_->{end} - $_->{start} foreach @$feat;
            $reghash->{$chr}=$feat;
        }
        return $reghash;
    }
}



