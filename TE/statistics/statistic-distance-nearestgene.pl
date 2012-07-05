use strict;
use warnings;
use Getopt::Long;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;


my $contigsonly=[qw/2L 2R 3L 3R 4 X/];
my $features=[qw/exon/];
my $contighash= {map {($_,1) } @$contigsonly};


my $teinsert_file;
my $annotation;
my $minfixratio;

GetOptions(
    "te-inserts=s"	    =>\$teinsert_file,
    "annotation=s"          =>\$annotation,
) or die "Wrong parameters";

die "insufficient parameters" unless -e $teinsert_file;

my $feath=Utility::load_features($annotation,$features,$contighash);
my $teinsraw=Utility::load_tes($teinsert_file);
my $tehash=Utility::createTEhash($teinsraw);


while(my($chr,$chrtes)=each(%$tehash))
{
    my $binchr=Utility::get_binary_chr_representation($feath->{$chr});
    foreach my $te(@$chrtes)
    {
        my $pos=int($te->{inspos});
        my $dowdist=0;
        my $usdist=0;
        
        my $i;
        for($i=$pos; $i>=0; $i--)
        {
            last if exists($binchr->[$i]);
            $dowdist++;
        }
        $dowdist=100000000 if $i==0; # lower chromosome boundary
        
        for($i=$pos;$i<=@$binchr; $i++)
        {
            last if exists($binchr->[$i]);
            $usdist++;
        }
        $usdist=100000000 if $i==@$binchr; # upper chromosome boundary


        my $dist=$dowdist < $usdist ? $dowdist : $usdist;
        print $te->{line}."\t$dist\n";
    }
    
    
    

    
}



exit;




{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    sub load_tes
    {
        my $file=shift;
        open my $ifh, "<", $file or die "Could not open input file";
        
        my $tes=[];
        while(my $l=<$ifh>)
        {
            chomp $l;
            #2L	10050000	10026010	FR	jockey	0.068859649	non-LTR	NA	119	10025851	10025950	0.05	80	4	76	0	10026070	10026138	0.087719298	57	5	52	0	4.03	euchr	60	522	1616.8	14127	intergenic	1	0.068859649	A	TRUE	In(2L)t	FALSE	NA	NA	NA	NA	NA	NA	5	5	0	0	0.166666667	2	2	2	2	1	1	0	0	0	0	0	0	0.2651515155	0.0798934935	0.482758621	27862	22164	353	35453	10000000	10100000	4.04	4.03	4.02	0.403	69136.476426799	54997.5186104218	E
            #2L	10050000	10033605.5	R	INE-1	0.482758621	TIR	FBti0064247	186	NA	NA	NA	NA	NA	NA	NA	10033804	10033887	0.482758621	29	14	15	0	4.03	euchr	40	152	253.9	7911	intron	1	0.482758621	A	TRUE	In(2L)t	FALSE	FBti0064247	2L	10033628	10033667	+	40	5	5	0	0	0.166666667	2	2	2	2	1	1	0	0	0	0	0	0	0.2651515155	0.0798934935	0.482758621	27862	22164	353	35453	10000000	10100000	4.04	4.03	4.02	0.403	69136.476426799	54997.5186104218	E
            my ($chr,undef,$pos)=split /\t/,$l;
            push @$tes,
            {
                chr=>$chr,
                inspos=>$pos,
                line=>$l
            };
        }
        return $tes;
    }
    
    sub get_length_from_binary
    {
        my $bin=shift;
        
        my $leng=0;
        foreach my $i (@$bin)
        {
            $leng++ if $i;
        }
        return $leng;
    }
    
    
    sub calc_intron
    {
        my $feath=shift;
        my $teh=shift;
        my $minfixratio=shift;
        
        my $feat_intr=$feath->{intron};
        my $feat_exon=$feath->{exon};
        
        my($intrl,$intrc,$intrfix)=(0,0,0);
        while(my($chr,$intron_feat)=each(%$feat_intr))
        {
            my $bin_chr=get_binary_chr_representation($intron_feat);
            
            my $exons=$feat_exon->{$chr};
            foreach my $ex (@$exons)
            {
                my $start=$ex->{start};
                my $end=$ex->{end};
                for my $i ($start..$end)
                {
                    $bin_chr->[$i]=0; 
                }
                
            }
            $intrl += get_length_from_binary($bin_chr);
            
            my $tes=$teh->{$chr};
            foreach my $te (@$tes)
            {
                my $hit=$bin_chr->[$te->{inspos}];
                $intrc++    if($hit);
                
                my $popfreq=$te->{popfreq}||0;
                $intrfix++  if($hit and $popfreq > $minfixratio);
            }
        }
        
        return ($intrl,$intrc,$intrfix);
    }
    
  
    
    
    sub get_binary_chr_representation
    {
        my $featurelist=shift;
        my $bin=[];
        my $chr_check="";
        foreach my $f (@$featurelist)
        {
            my $chr=$f->{chr};
            $chr_check=$chr unless $chr_check;
            die "chromosome fucked" unless $chr_check eq $chr;
            my $start=$f->{start};
            my $end=$f->{end};
            for my $i ($start..$end)
            {
                $bin->[$i]=1;
            }
        }
        return $bin;
    }
    
    
    sub createTEhash
    {
        my $tes=shift;
        
        my $teh;
        foreach my $te (@$tes)
        {
            # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
            # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
            # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
            my $chr=$te->{chr};
            $teh->{$chr}=[] unless(exists($teh->{$chr}));
            push @{$teh->{$chr}},$te;
        }
        return $teh;
    }
    
    
    
    sub load_features
    {
        my $file        =   shift;
        my $features    =   shift;
        my $contighash  =   shift;
        my $featurehash =   {map {($_,1) } @$features};
        
        open my $ifh, "<", $file or die "could not open input file";
        
        
        my $toret={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            next if $l=~/^#/;
            my $p = parseFlybase($l);
            my $cat=$p->{cat};
            my $chr=$p->{chr};
            next unless(exists($featurehash->{$cat}));
            next unless(exists($contighash->{$chr}));
            $toret->{$chr} = [] unless(exists($toret->{$chr}));
            push @{$toret->{$chr}},$p;
        }
        return $toret;
    }
    
    sub parseFlybase
    {
        my $l=shift;
        chomp $l;
        #2L	FlyBase	chromosome_band	-204333	22221	.	+	.	ID=band-21A_chromosome_band;Name=band-21A
        #2L	FlyBase	chromosome_band	-204333	-153714	.	+	.	ID=band-21A1_chromosome_band;Name=band-21A1
            
        my($chr,undef,$cat,$start,$end,undef,$strand,undef,$misc)=split /\t/,$l;
        my $e =
        {
              chr=>$chr,
              cat=>$cat,
              start=>$start,
              end=>$end,
              strand=>$strand,
              misc=>$misc
        };
        return $e;
    }
    
 
    
    
    

    
    

    
 
}





