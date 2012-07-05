use strict;
use warnings;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use List::Util qw[min max];
use TEInsertUtility;
use Recombination;
use Getopt::Long;


my $teinsertfile="";
my $outputfile="";

my $telengthfile="";
my $recratefile="";
my $recratethreshold=1;

my $annotationfile="";


GetOptions(
    "teinsert=s"	    =>\$teinsertfile,
    "output=s"              =>\$outputfile,
    "annotation=s"          =>\$annotationfile,
    "teleng-file=s"         =>\$telengthfile,
    "recombination-file=s"        =>\$recratefile,
    "recombination-threshold=s"   =>\$recratethreshold,
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

my $recrate_annotator=Utility::get_recrate_annotator($recratefile);
my $chrregion_annotator=Utility::get_chrregion_annotator($recratefile,$recratethreshold);
my $teleng_annotator=Utility::get_telength_annotator($telengthfile);
my $annotation_annotator=Utility::get_annotation_annotator($annotationfile);


open my $ofh,">",$outputfile or die "Could not open output file";
my $polyins=load_te_inserts($teinsertfile);
foreach my $ins (@$polyins)
{
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap

    my $teid = $ins->{teid};
    my $chr  =  $ins->{chr};
    my $pos  = int($ins->{inspos});

    
    my $recrate=$recrate_annotator->($chr,$pos);
    my $chrregion=$chrregion_annotator->($chr,$pos);
    my $famleng=$teleng_annotator->($teid);
    my $annot=$annotation_annotator->($chr,$pos);
    
    if(not defined($recrate) or not defined($chrregion) or not defined($famleng) or not defined($annot))
    {
        my $bla=0;
    }
    print $ofh $ins->tostring()."\t$recrate\t$chrregion\t$famleng\t$annot\n"; 


}

exit;




{
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use List::Util qw[min max];
    use Recombination;
    package Utility;
    
    
    sub get_recrate_annotator
    {
        my $recratefile=shift;
        my $rr=Recombination::load_recombination_rate($recratefile);
        my $windowsize=100000;
        
        return sub
        {
            my $chr=shift;
            my $pos=shift;
            my $index=int($pos/$windowsize);
            my $reckey=((2*$index+1)*$windowsize)/2;
            my $rec=$rr->{$chr}{$reckey};
            $rec||=0;
            return $rec;
            
        }
    }
    
    sub get_telength_annotator
    {
        my $telengthfile=shift;
        my $famleng={};
        open my $ifh, "<",$telengthfile or die "Could not open input file";
        while(my $l=<$ifh>)
        {
            chomp $l;
            my($fam,$leng)=split /\t/,$l;
            $famleng->{$fam}=$leng;
        }
        
        return sub{
            my $queryfam=shift;
            die "no entry for family $queryfam" unless(exists($famleng->{$queryfam}));
            return $famleng->{$queryfam};
        }
    }
    
    
    sub get_chrregion_annotator
    {
        my $recombinationfile=shift;
        my $recthreshold=shift;
        my $chrinfo={
                    X       =>{leng=>22422827, cent_pos=>"end"},
                    '2L'    =>{leng=>23011544, cent_pos=>"end"},
                    '2R'    =>{leng=>21146708, cent_pos=>"start"},
                    '3L'    =>{leng=>24543557, cent_pos=>"end"},
                    '3R'    =>{leng=>27905053, cent_pos=>"start"}
                    };
        
        my $regh=_define_regions($chrinfo,$recombinationfile,$recthreshold);
        
        return sub{
            my $chr=shift;
            my $pos=shift;
            my $cregs=$regh->{$chr};
            foreach my $ar (@$cregs)
            {
                if($pos>$ar->{start} and $pos<=$ar->{end})
                {
                    my $feature=$ar->{feat};
                    unless($feature)
                    {
                        die "Could not find feature for $chr $pos";    
                    }
                    return $feature;
                }
            }
            return "na";   
        }
    }
    
    sub get_annotation_annotator
    {
        my $annotationfile=shift;
        my $annotationlist=__load_annotation($annotationfile);
        return sub
        {
            my $chr=shift;
            my $pos=shift;
            my $feat=_search_feature($annotationlist,$chr,$pos);
            die "No feature $feat for $chr $pos" unless $feat;
            return $feat;
        }
    }
    
    sub __load_annotation{
        my $annotationfile=shift;
        my $contighash={"X"=>1,"2L"=>1,"2R"=>1,"3L"=>1,"3R"=>1,"4"=>1};
        my $catcoder={chromosome_arm=>"n", CDS=>"c", exon=>"e",five_prime_UTR=>"5",intron=>"i",pseudogene=>"p",three_prime_UTR=>"3",regulatory_region=>"r"};
        open my $ifh,"<",$annotationfile or die "Could not open annotation file";
        
        # load annotation
        my $annotationlist={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            my $p=__parseFlybase($l);
            my $chr=$p->{chr};
            my $cat=$p->{cat};
            next unless(exists($contighash->{$chr}));
            next unless(exists($catcoder->{$cat}));
            $annotationlist->{$chr}{$cat}||=[];
            push @{$annotationlist->{$chr}{$cat}},$p;
        }

        return $annotationlist;
    }
    
    
    sub _search_feature
    {
        my $featurelist=shift;
        my $chr=shift;
        my $pos=shift;
        
        my $feat="";
        $feat=__search_featurelist($featurelist->{$chr}{CDS},$pos);
        return "CDS" if $feat;
        $feat=__search_featurelist($featurelist->{$chr}{three_prime_UTR},$pos);
        return "3-UTR" if $feat;
        $feat=__search_featurelist($featurelist->{$chr}{five_prime_UTR},$pos);
        return "5-UTR" if $feat;
        $feat=__search_featurelist($featurelist->{$chr}{exon},$pos);
        return "exon" if $feat;
        $feat=__search_featurelist($featurelist->{$chr}{pseudogene},$pos);
        return "pseudogene" if $feat;
        $feat=__search_featurelist($featurelist->{$chr}{regulatory_region},$pos);
        return "regulatory" if $feat;
        $feat=__search_featurelist($featurelist->{$chr}{intron},$pos);
        return "intron" if $feat;
        $feat=__search_featurelist($featurelist->{$chr}{chromosome_arm},$pos);
        return "intergenic" if $feat;
        return "na";
    }
    
    sub __search_featurelist
    {
        my $featurelist=shift;
        my $pos=shift;
        
        die "no position" unless $pos;
        return undef unless $featurelist;
        foreach my $feat (@$featurelist)
        {
            my $start=$feat->{start};
            my $end=$feat->{end};
            if($pos>=$start and $pos <=$end)
            {
                return $feat->{cat};
            }
        }
        return undef;
 


    }
    
    sub __parseFlybase
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
    
    
    sub __load_recombination_rate
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
    
    sub _define_regions
    {
        my $chrinfo=shift;
        my $recombinationfile=shift;
        my $recthreshold=shift;
        
        my $rh=__load_recombination_rate($recombinationfile);
        
        
        
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