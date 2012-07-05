use strict;
use warnings;
use Getopt::Long;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;


my $contigsonly=[qw/2L 2R 3L 3R 4 X/];
my $features=[qw/exon CDS/];
my $contighash= {map {($_,1) } @$contigsonly};



my $annotation;
my $stepsize;

GetOptions(

    "annotation=s"          =>\$annotation,
    "step-size=i"           =>\$stepsize

) or die "Wrong parameters";

die "provide a annotation" unless -e $annotation;
die "provide a step-size" unless $stepsize;

my $feath=Utility::load_features($annotation,$features,$contighash);



while(my($chr,$temp)=each(%$feath))
{

    foreach my $cat (@$features)
    {
        my $start=1;
        my $end=$stepsize;
        my $featurelist=$temp->{$cat};
        my $binfeat=Utility::get_binary_chr_representation($featurelist);
        my $length=@$binfeat;
        
        while($start<$length)
        {
            my $count=0;
            for my $i($start..($end-1))
            {
                $count++ if($binfeat->[$i]);
            }
            
            my $pos=int(($start+$end)/2);
            
            print "$cat\t$chr\t$pos\t$count\n";
            
            $start+=$stepsize;
            $end+=$stepsize;
        }

    }
    
    while(my($chr,$featurelist)=each(%$temp))
    {


    }
}


exit;




{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
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
            $toret->{$cat}{$chr} = [] unless(exists($toret->{$chr}{$cat}));
            push @{$toret->{$chr}{$cat}},$p;
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





