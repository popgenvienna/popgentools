use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use List::Util qw[min max];
use Getopt::Long;
use Pod::Usage;

my $generations =0;
my $popSize     =0;
my $help        =0;
my $synchronizedfile="";
my $mincount=2;
my $mincoverage=4;
my $maxcoverage=400;
my $printonlysnps=0;

GetOptions(
    "min-count=i"               =>\$mincount,
    "min-coverage=i"            =>\$mincoverage,
    "max-coverage=i"            =>\$maxcoverage,
    "synchronized=s"            =>\$synchronizedfile,
    "generations=s"	        =>\$generations,
    "pop-size=i"                =>\$popSize,
    "only-snps"                 =>\$printonlysnps,
    "help"                      =>\$help            
    ) or die "unknown option\n";

pod2usage(-msg=>"provide synchronized file") unless $synchronizedfile;
pod2usage(-msg=>"provide --generations") unless $generations;
pod2usage(-msg=>"provide --pop-size") unless $popSize;
pod2usage(-verbose=>2) if $help;


my $gen=[split /:/,$generations];


open my $ifh, "<",$synchronizedfile or die "could not open file";


while(my $l =<$ifh>)
{
    chomp $l;
    my $p=Utility::parseSync($l);
    if($p->{mincov}<$mincoverage or $p->{maxcov}>$maxcoverage or $p->{sec_al_count} < $mincount or $p->{del} > 0)
    {
        Utility::printsync($p) unless $printonlysnps;
        next;
    }
    
    # let the mutations begin
    my $base=$p->{base};
    my $derived=$p->{derived};
    my ($cgen,$cbase,$cderived) = (scalar(@$gen),scalar(@$base),scalar(@$derived));
    die "the vector specifying the number of generations has a different size as the number of base populations" unless $cgen == $cbase;
    die "the number of base populations does not agree with the number of derived populations" unless $cbase== $cderived;
    
    for my $i(0..($cbase-1))
    {
        my $act_ba=$base->[$i];
        my $act_de=$derived->[$i];
        my $act_gen=$gen->[$i];
        my $target_cov=$act_de->{eucov};
        
        my $major=$act_ba->{major};
        my $minor=$act_ba->{minor};
        my $majorcount=$act_ba->{$major};
        my $minorcount=$act_ba->{eucov}-$majorcount;
        $act_ba->{A}=0;$act_ba->{T}=0;$act_ba->{C}=0;$act_ba->{G}=0;
        $act_ba->{$major}=$majorcount;$act_ba->{$minor}=$minorcount;
        
        my($der_minc,$der_majc) = Utility::get_derived_counts($majorcount,$minorcount,$target_cov,$popSize,$act_gen);
    

        $act_de->{A}=0; $act_de->{T}=0;$act_de->{C}=0;$act_de->{G}=0;$act_de->{N}=0; $act_de->{del}=0;
        $act_de->{$major}=$der_majc; $act_de->{$minor}=$der_minc;
    }
    Utility::printsync($p);
}

exit;

{
    package Utility;
    use strict;
    use warnings;
    

    
    sub parseSync
    {
    #2L	79	G	0:0:0:15:0:0	0:0:0:38:0:0
    #2L	80	A	12:0:0:0:0:0	38:0:0:0:0:0
    #2L	81	A	14:0:0:0:0:0	43:0:0:0:0:0
    #2L	82	A	14:0:0:0:0:0	42:0:0:0:0:0
        my $line=shift;
        chomp $line;
        my @ar=split /\s+/,$line;
        my $chr =shift @ar;
        my $pos= shift @ar;
        my $refc=shift @ar;
        
        my @pops=map {_parseSyncEntry($_)} @ar;
        my $popcount=@pops;
        die "number of populations in the sync file has to be an even number" unless ($popcount%2==0);
        my $halfpop=int($popcount/2);
        
        
        ################# get min max cov##############
        my $mincov=1000000000;
        my $maxcov=0;
        foreach my $p (@pops)
        {
            $maxcov=$p->{eucov} if($p->{eucov} > $maxcov);
            $mincov=$p->{eucov} if($p->{eucov} < $mincov);
        }
        
        ######### separate base population from derived populations
        my @base=();
        for my $i (0..($halfpop-1))
        {
            my $e= shift @pops;
            push @base,$e;
        }
        
        
        for my $i (0..(scalar(@pops)-1))
        {
            my $ba=$base[$i];
            my $der=$pops[$i];
            my $bacov=$ba->{eucov};
            my $dercov=$der->{eucov};
            next if $bacov==0;
            next if $dercov==0;
            my $scaling=$dercov/$bacov;
            $der->{A}=int($ba->{A}*$scaling);
            $der->{T}=int($ba->{T}*$scaling);
            $der->{C}=int($ba->{C}*$scaling);
            $der->{G}=int($ba->{G}*$scaling);
            $der->{N}=int($ba->{N}*$scaling);
            $der->{del}=int($ba->{del}*$scaling);

        }
       
        
        my($ac,$tc,$cc,$gc,$nc,$dc)=(0,0,0,0,0,0);
        foreach my $e (@base)
        {
            $ac+=$e->{A};
            $tc+=$e->{T};
            $cc+=$e->{C};
            $gc+=$e->{G};
            $nc+=$e->{N};
            $dc+=$e->{del};
        }
        
        my $secalcount=0;
        my @tar=($ac,$tc,$cc,$gc);
        @tar = sort {$b<=>$a} @tar;
        $secalcount=$tar[1];
        return  {
            chr=>$chr,
            pos=>$pos,
            refc=>$refc,
            base=>\@base,
            derived=>\@pops,
            mincov=>$mincov,
            maxcov=>$maxcov,    
            del=>$dc,
            sec_al_count=>$secalcount
        };
    }
    
    sub printsync
    {
        my $p=shift;
        my $head="$p->{chr}\t$p->{pos}\t$p->{refc}";
        my $base=$p->{base};
        my $derived=$p->{derived};
        
        my @ar=(@$base,@$derived);
        my @temp=();
        foreach my $e (@ar)
        {
            my $te="$e->{A}:$e->{T}:$e->{C}:$e->{G}:$e->{N}:$e->{del}";
            push @temp,$te;
        }
        my $str=join("\t",@temp);
        print "$head\t$str\n";
    }
    
    sub _parseSyncEntry
    {
        my $l=shift;
        
        return {A=>0,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>0,major=>"N",minor=>"N"} if $l eq "-";
        my($ac,$tc,$cc,$gc,$nc,$dc)=split /:/,$l;
        
        my $eucov=$ac+$tc+$gc+$cc;
        my @als=({a=>"A",c=>$ac},{a=>"T",c=>$tc},{a=>"C",c=>$cc},{a=>"G",c=>$gc});
        @als=sort {$b->{c}<=>$a->{c}} @als;
        my $maj=$als[0]->{a};
        my $minor=$als[1]->{a};
        
        return {
            A   =>$ac,
            T   =>$tc,
            C   =>$cc,
            G   =>$gc,
            del =>$dc,
            N   =>$nc,
            eucov=>$eucov,
            major=>$maj,
            minor=>$minor
        };
    }
    
    
    # MINOR ALLELE = A
    # MAJOR ALLELE = C
    
    sub get_derived_counts
    {
        my $majorcount =shift;
        my $minorcount =shift;
        my $target_cov=shift;
        my $popSize=shift;
        my $generations=shift;
        
        my $ful_cov=$majorcount+$minorcount;
        
        my $targetfreq=$minorcount/$ful_cov;
        # minor allele=A; major allele = C
        my $basepop = get_basepopulation($targetfreq,$popSize);
        my $drifted = drift_population($basepop,$generations);
        my $covered= random_draw ($drifted,$target_cov);
        my $af=measure_allelefreq($covered);
        return ($af->{A},$af->{C});
    }

    
    sub measure_allelefreq
    {
        my $pop=shift;
        my $af={};
        my $cov=0;
        foreach my $a (@$pop)
        {
            $cov++;
            $af->{$a}++;
        }
        $af->{A}||=0;
        $af->{C}||=0;
        $af->{cov}=$cov;
        return $af;
    }
    
    sub random_draw
    {
        # draw randomly some alleles from the total; this is simulating the coverage from pooled data
        my $pop=shift;
        my $avcov=shift;
        
        my $toret=[];
        
        for (1..$avcov)
        {
            my $a=$pop->[int(rand(scalar(@$pop)))];
            push @$toret,$a;
        }
        return $toret;
    }
    


    sub get_threshold
    {
        # calculate the minor allele frequency: aka threshold
        # minor allele = A
        my $pop=shift;
        my $a_count=0;
        my $cov=0;
        foreach my $a (@$pop)
        {
            $cov++;
            $a_count++ if $a eq "A";
        }
        my $af=$a_count/$cov;
        return $af;
    }


    
    sub get_basepopulation
    {
        my $targetfreq  =shift;
        my $popSize     =shift;
        my $al_a = int($popSize*2*$targetfreq); # minor allele = A
        my $al_b = 2 * $popSize - $al_a;        # major allele = C
        
        my $pop= [(split //, "A" x $al_a ), (split //,"C" x $al_b) ];
        return $pop;
    }
    
    
    
    sub drift_population
    {
        my $basepop=shift;
        my $generation=shift;
        
        my $threshold=get_threshold($basepop);
        my $drift_pop;
        for (1..$generation)
        {
            # drift the population
            $drift_pop=_get_next_generation(scalar(@$basepop),$threshold);
            # and recalculate the threshold for the next generation
            $threshold=get_threshold($drift_pop);
        }
        return $drift_pop;
    }
    
    sub _get_next_generation
    {
        my $pop_size=shift;
        my $threshold=shift;
        
        my $driftp=[];
        for (1..$pop_size)
        {
            my $val=rand;
            # now this is important
            # the minor allele is A;
            # the major allele is C
            # when calculating the threshold the A are counted
            # so if the random variable is smaller or equal to the threshold it should be the minor allele A
            # if larger it should be C
            my $a=$val>$threshold?"C":"A";
            push @$driftp,$a;
        }
        return $driftp;
    }
}

    #"min-count=i"               =>\$mincount,
    #"min-coverage=i"            =>\$mincoverage,
    #"max-coverage=i"            =>\$maxcoverage,
    #"synchronized=s"            =>\$synchronizedfile,
    #"generations=s"	        =>\$generations,
    #"pop-size=i"                =>\$popSize,
    #"only-snps"                 =>\$printonlysnps,
    #"help"                      =>\$help        

=head1 NAME

simulate-drift-multibasepop.pl - simulates drift and creates a new synchronized file containing the drifted populations

=head1 SYNOPSIS
    
perl simulate-drift-multibasepop.pl --synchronized input.syn --generations 23:35 --pop-size 500 --min-count 4 --min-coverage 100 --max-coverage 400 --only-snps

=cut

