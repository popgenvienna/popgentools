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
use Getopt::Long;

my $featfile=shift;
my $nfile=shift;
my $windowsize=shift;

my $feath = Utility::load_features($featfile);
my $nh=Utility::load_ngtf($nfile);

while( my($chr,$featurelist) = each(%$feath))
{
    my $start=1;
    my $end=$windowsize;
    
    my $nlist=$nh->{$chr}||[];
    $featurelist->{te}=$nlist;
    
    my $binchr=Utility::get_intergenic_binary($featurelist);
    my $length=@$binchr;
    while($start<$length)
    {
        my $ncount=0;
        for my $i (($start)..$end)
        {
            $ncount++ if $binchr->[$i];
        }
        

        my $mid=int(($start+$end)/2);
        print "$chr\t$mid\t$ncount\n";
        
        $start+=$windowsize;
        $end+=$windowsize;
    }
    
}




{
    package Utility;
    use strict;
    use warnings;
    
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
    
    sub load_ngtf
    {
        my $file        =   shift; 
        my $contighash  =   {"2L"=>1,"X"=>1,"2R"=>1,"3L"=>1,"3R"=>1,"4"=>1};
        my $featurehash =   {polyN => 1};
        
        open my $ifh, "<", $file or die "could not open input file";
        # YHet    polyNsearch     polyN   5917    5991    .       +       .       gene_id "poly_N_9";
        # YHet    polyNsearch     polyN   6002    13762   .       +       .       gene_id "poly_N_10";

        
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
    
    sub get_intergenic_binary
    {
        my $feath=shift;
        
        my $chromosome_arm=$feath->{chromosome_arm};
        
    
        my $bin_chr=get_binary_chr_representation($feath->{chromosome_arm});
        while(my($cat,$list)=each(%$feath))
        {
            next if $cat eq "chromosome_arm";
            foreach my $feat (@$list)
            {
                my $start=$feat->{start};
                my $end=$feat->{end};
                for my $i ($start..$end)
                {
                    $bin_chr->[$i]=0;  
                }
            }
                    
        }
        return $bin_chr;
    }
    
    
    
    
    sub load_features
    {
        my $file        =   shift; 
        my $contighash  =   {"2L"=>1,"X"=>1,"2R"=>1,"3L"=>1,"3R"=>1,"4"=>1};
        my $featurehash =   {intron=>1, exon=>1, CDS=>1, five_prime_UTR=>1, three_prime_UTR=>1, ncRNA=>1, regulatory_region=>1, chromosome_arm=>1};
        
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
            $toret->{$chr}{$cat} = [] unless(exists($toret->{$chr}{$cat}));
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