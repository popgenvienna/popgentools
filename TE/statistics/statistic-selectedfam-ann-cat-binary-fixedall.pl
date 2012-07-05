use strict;
use warnings;
use Getopt::Long;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEInsertUtility;


my $contigsonly=[qw/2L 2R 3L 3R 4 X/];
my $features=[qw/intron exon CDS five_prime_UTR three_prime_UTR ncRNA regulatory_region chromosome_arm/];
my $contighash= {map {($_,1) } @$contigsonly};


my $teinsert_file;
my $annotation;
my $minfixratio;
my $selectedfamilyfile;

GetOptions(
    "te-inserts=s"	    =>\$teinsert_file,
    "annotation=s"          =>\$annotation,
    "min-fixratio=f"        =>\$minfixratio,
    "selected-family=s"     =>\$selectedfamilyfile
) or die "Wrong parameters";

die "insufficient parameters" unless -e $teinsert_file;
die "wrong minfixratio" unless $minfixratio;

my $families=Utility::load_selected_family($selectedfamilyfile);
my $feath=Utility::load_features($annotation,$features,$contighash);
my $teinsall=load_te_inserts($teinsert_file);
my $teinsraw=Utility::filter_selectedfamiles($teinsall,$families);
my $teh=Utility::createTEhash($teinsraw);

my $annLen={};
my $annCount={};
my $fixCount={};

my($intgenleng,$intgenCount,$intgenFix)=Utility::calc_intergenic($feath,$teh,$minfixratio);
my($intronleng,$intronCount,$intronFix)=Utility::calc_intron($feath,$teh,$minfixratio);

while(my($cat,$temp)=each(%$feath))
{
    while(my($chr,$featurelist)=each(%$temp))
    {
        my $telist=$teh->{$chr};
        my $bin_chr=Utility::get_binary_chr_representation($featurelist);
        my $length=Utility::get_length_from_binary($bin_chr);
        $annLen->{$cat}+=$length;
        foreach my $te (@$telist)
        {
            my $hit=$bin_chr->[$te->{inspos}];
            $annCount->{$cat}++ if $hit;
            my $popfreq=$te->{popfreq};
            $fixCount->{$cat}++ if($hit and $popfreq> $minfixratio);
        }
    }
}

while(my($cat,$leng)=each(%$annLen))
{
    my $allcount=$annCount->{$cat} ||0;
    my $fixcount=$fixCount->{$cat} ||0;
    my $alldensity=sprintf("%.1f",$allcount*1000000/$leng);
    my $fixdensity=sprintf("%.1f",$fixcount*1000000/$leng);
        
    
    print "$cat\t$leng\t$allcount\t$fixcount\t$alldensity\t$fixdensity\n";
}

my $intgdens=sprintf("%.1f",$intgenCount*1000000/$intgenleng);
my $fixdens=sprintf("%.1f",$intgenFix*1000000/$intgenleng);
print "intergenic\t$intgenleng\t$intgenCount\t$intgenFix\t$intgdens\t$fixdens\n";

my $introndens=sprintf("%.1f",$intronCount*1000000/$intgenleng);
my $intronfixdens=sprintf("%.1f",$intronFix*1000000/$intgenleng);
print "intron(excl-exon)\t$intronleng\t$intronCount\t$intronFix\t$introndens\t$intronfixdens\n";
exit;




{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    
    sub filter_selectedfamiles
    {
        my $teins=shift;
        my $filterfam=shift;
        my $toret=[];
        foreach my $te (@$teins)
        {
            my $fam=$te->{teid};
            next unless(exists($filterfam->{$fam}));
            push @$toret,$te;
        }
        return $toret;
    }
    
    sub load_selected_family
    {
        my $file=shift;
        open my $ifh, "<", $file or die "Could not open input file";
        
        my $families={};
        while(my $line=<$ifh>)
        {
            chomp $line;
            $families->{$line}=1
        }
        return $families;
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
    
    sub calc_intergenic
    {
        my $feath=shift;
        my $teh=shift;
        my $minfixratio=shift;
        
        my $chromosome_arm=$feath->{chromosome_arm};
        
        my($intgl,$intgc,$intgfix)=(0,0,0);
        while(my($chr,$chrarmfeat)=each(%$chromosome_arm))
        {
            my $bin_chr=get_binary_chr_representation($chrarmfeat);
            while(my($cat,$temp)=each(%$feath))
            {
                next if $cat eq "chromosome_arm";
                my $tempfeat=$temp->{$chr};
                foreach my $feat (@$tempfeat)
                {
                    my $start=$feat->{start};
                    my $end=$feat->{end};
                    for my $i ($start..$end)
                    {
                        $bin_chr->[$i]=0;  
                    }
                    
                }
            }
            $intgl += get_length_from_binary($bin_chr);
            
            
            my $tes=$teh->{$chr};
            foreach my $te (@$tes)
            {
                my $hit=$bin_chr->[$te->{inspos}];
                $intgc++    if($hit);    
                my $popfreq=$te->{popfreq}||0;
                $intgfix++  if($hit and $popfreq > $minfixratio);
            }
        }
        
        return ($intgl,$intgc,$intgfix);
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
            $toret->{$cat}{$chr} = [] unless(exists($toret->{$cat}{$chr}));
            push @{$toret->{$cat}{$chr}},$p;
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





