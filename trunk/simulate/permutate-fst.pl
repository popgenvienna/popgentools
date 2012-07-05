#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; 
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use Synchronized;



my $input;
my $output;
my $help=0;
my $test=0;
my $popsize=0;
my $classcount=10;
my $permutcount=10000;
my $mincount=0;
my $mincoverage=0;
my $maxcoverage=0;
my $rawchromosomes="";


# --input /Volumes/Volume_4/analysis/pablo-permut/test.syn --classes 100 --permutations 10000 --popsize 500
GetOptions(
    "input=s"	        =>\$input,
    "output=s"          =>\$output,
    "classes=i"         =>\$classcount,
    "permutations=i"    =>\$permutcount,
    "min-count=i"       =>\$mincount,
    "min-coverage=i"    =>\$mincoverage,
    "max-coverage=i"    =>\$maxcoverage,
    "chromosomes=s"     =>\$rawchromosomes,
    "popsize=i"         =>\$popsize,
    "test"              =>\$test,
    "help"	        =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);
pod2usage(-verbose=>2) if $help;
pod2usage() unless -e $input;
pod2usage() unless $classcount>1;
pod2usage() unless $output;
pod2usage() unless $mincount;
pod2usage() unless $popsize;
pod2usage() unless $mincoverage;
pod2usage() unless $maxcoverage;

my @chrar=split/\s+/,$rawchromosomes;
push @chrar,$rawchromosomes unless $rawchromosomes;

my $chrhash= {map {($_,1)} @chrar};
my $sp=get_sumsnp_synparser($mincount,$mincoverage,$maxcoverage);

open my $ifh,"<", $input or die "Could not open input file";

# all SNPs are stored in an array, so they can be randomly chosen
my $s_all=[];
my $s_class=[]; # and all SNPs are categorized into classes;

while(my $line=<$ifh>)
{
    my $p=$sp->($line);
    # throw all non-pure SNP away;
    next unless $p->{ispuresnp};
    next unless(exists($chrhash->{$p->{chr}}));
    
    my $rc=$p->{refchar};
    my $samples=$p->{samples};
    my $base=$p->{samples}[0];
    
    # count all reference bases in the populations and discard if smaller than mincount
    my @refcounts=map {$_->{$rc}} @$samples;
    my $totrefcount=0;
    $totrefcount+=$_ foreach @refcounts;
    next unless $totrefcount>=$mincount;

    # calculate the frequency of the reference base in the base population and calculate the proper category
    my $basereffreq = $base->{$rc} / $base->{eucov};
    my $cat=int($classcount*$basereffreq);
    $p->{cat}=$cat;
    $p->{reff}=$basereffreq;
    
    # add the SNP to the collection
    push @$s_all,$p;
    $s_class->[$cat]=[] unless $s_class->[$cat];
    push @{$s_class->[$cat]},$p;
}

# Print summary statistics

my $count_allsnps=@$s_all;
my $count_classes =@$s_class;
print "Using $count_allsnps SNPs; Distributed to $classcount classes\n";
for my $i (0..$classcount)
{
    my $temp=$s_class->[$i];
    if($temp)
    {
        $temp=@$temp;
    }
    else
    {
        $temp=0;
    }
    print "Class: $i - Count: $temp\n";
}

open my $ofh, ">", $output or die "Could not open output file";

for my $i (1..$permutcount)
{
    # pick a random snp
    my $temp=int(rand(@$s_all));
    my $e_base=$s_all->[$temp];

    # get all SNPs from the same category;
    my $cat=$e_base->{cat};
    my $alcat=$s_class->[$cat];
    my $count_incat=@$alcat;
    
    #permutate in the same class
    my $e_p2=Utility::get_random_subpop($alcat);
    my $e_p3=Utility::get_random_subpop($alcat);
    my $e_p4=Utility::get_random_subpop($alcat);
    
    # normalize the sample
    my $base=Utility::normalize_sample($e_base->{samples}[0],$e_base->{refchar});
    my $p2=Utility::normalize_sample($e_p2->{samples}[1],$e_p2->{refchar});
    my $p3=Utility::normalize_sample($e_p3->{samples}[2],$e_p3->{refchar});
    my $p4=Utility::normalize_sample($e_p4->{samples}[3],$e_p4->{refchar});
    
    
    my($pitot12,$fst12) = Utility::calculateFst($base,$p2,$popsize);
    my($pitot13,$fst13) = Utility::calculateFst($base,$p3,$popsize);
    my($pitot14,$fst14) = Utility::calculateFst($base,$p4,$popsize);
    my($pitot23,$fst23) = Utility::calculateFst($p2,$p3,$popsize);
    my($pitot34,$fst34) = Utility::calculateFst($p3,$p4,$popsize);
    my($pitot24,$fst24) = Utility::calculateFst($p2,$p4,$popsize);
    
    
    
    print $ofh "$i\t1:2\t$fst12\t$base->{Ac}\t$base->{Cc}\t$p2->{Ac}\t$p2->{Cc}\t$pitot12\t$cat\t$count_incat\n";
    print $ofh "$i\t1:3\t$fst13\t$base->{Ac}\t$base->{Cc}\t$p3->{Ac}\t$p3->{Cc}\t$pitot13\t$cat\t$count_incat\n";
    print $ofh "$i\t1:4\t$fst14\t$base->{Ac}\t$base->{Cc}\t$p4->{Ac}\t$p4->{Cc}\t$pitot14\t$cat\t$count_incat\n";
    print $ofh "$i\t2:3\t$fst23\t$p2->{Ac}\t$p2->{Cc}\t$p3->{Ac}\t$p3->{Cc}\t$pitot23\t$cat\t$count_incat\n";
    print $ofh "$i\t3:4\t$fst34\t$p3->{Ac}\t$p3->{Cc}\t$p4->{Ac}\t$p4->{Cc}\t$pitot34\t$cat\t$count_incat\n";
    print $ofh "$i\t2:4\t$fst24\t$p2->{Ac}\t$p2->{Cc}\t$p4->{Ac}\t$p4->{Cc}\t$pitot24\t$cat\t$count_incat\n";
    # i 1:2 Fst A1 C1 A2 C2   
    
}

exit;


{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    
    sub get_random_subpop
    {
        # provided all SNPs from a certain category;
        # this function selects a random SNP and returns the sample with the provided ID
        my $cat=shift;
        my $sampleid=shift;
        
        my $count=@$cat;
        my $rnum=int(rand($count));
        
        my $rsamp=$cat->[$rnum];
        return $rsamp;
    }
    
    sub normalize_sample
    {
        # create a fake SNP based on the data of the true SNP; A is always the reference allele and C always the alternative allele
        my $sample=shift;
        my $rc=shift;
        my $cov=$sample->{eucov};
        my $refcount=$sample->{$rc};
        
        # count of alternative allele
        my $altcount=$cov-$refcount;
        my $reffreq=$refcount/$cov;
        return
        {
            A=>$reffreq,
            Ac=>$refcount,
            C=>(1-$reffreq),
            Cc=>$altcount,
            cov=>$cov
        };
    }
    
    sub calculateFst
    {
        # calculte the Fst; provided are the allele counts for two populations
        my $pop1=shift;
        my $pop2=shift;
        my $popSize=shift;
        
        my $totpop=
        {
            A=>($pop1->{A}+$pop2->{A})/2, # dangerous
            C=>($pop1->{C}+$pop2->{C})/2,
            cov=>min($pop1->{cov},$pop2->{cov})
        };
        
        my $pip1=_pi($pop1,$popSize);
        my $pip2=_pi($pop2,$popSize);
        my $pitot=_pi($totpop,$popSize);
        
        my $av=($pip1+$pip2)/2;
        
        my $fst=0;
        $fst=(($pitot-$av)/$pitot) if $pitot;
        return ($pitot,$fst);
    }
    
    sub _pi
    {
        my $snp = shift;
        my $popSize=shift;
        my $M=$snp->{cov};
        my $pi_snp=1;
        $pi_snp-=(($snp->{A}))**2;
        $pi_snp-=(($snp->{C}))**2;
    
        #this term is indepent of the correction
        $pi_snp*=($M/($M-1));
        $pi_snp*=($popSize/($popSize-1));
        return $pi_snp;
    }
}


# Problems:
#1. mincoverage, for Fst we take the min coverage, the coverage between several populations may actually be correlated due to GC bias
# we are ignoring this thing totally by permutating;
#2. I have to use two fake alleles as SNPs with different major and minor alleles will enter the same categories; if not we may get unreasonable high Fst
# I call the reference allele A and the alternative allele C;
#3. A is always the frequency of the reference allele and C always 1 - ref_allele_freq -> eliminating three allelic states
#4. In addition to the SNP criteria, the reference allele has to fullfill the SNP criteria as well (mincount) -> otherwise unrealistic low fsts

#2L	5762	T	0:69:29:0:0:0	0:20:7:0:0:0	0:46:11:1:0:0	0:26:14:0:0:0	0:51:16:0:0:0
#2L	5776	C	0:0:96:3:0:0	0:0:19:0:0:0	0:0:45:0:0:0	0:0:27:0:0:0	0:0:71:0:0:0
#2L	5808	A	96:1:0:0:0:0	22:1:0:0:0:0	65:0:0:0:0:0	44:0:0:0:0:0	69:0:0:0:0:0
#2L	5813	G	0:4:0:94:0:0	0:3:0:12:0:0	0:1:0:41:0:0	0:5:0:37:0:0	0:1:0:65:0:0

#GetOptions(
#    "input=s"	        =>\$input,
#    "output=s"          =>\$output,
#    "classes=i"         =>\$classcount,
#    "permutations=i"    =>\$permutcount,
#    "min-count=i"       =>\$mincount,
#    "min-coverage=i"    =>\$mincoverage,
#    "max-coverage=i"    =>\$maxcoverage,
#    "popsize=i"         =>\$popsize,
#    "test"              =>\$test,
#    "help"	        =>\$help

=head1 NAME

permutate-fst.pl - creates permutation tests for Martin and Pablo; 

=head1 SYNOPSIS

 perl permutate-fst.pl --popsize 500 --input exp.sync --output permutated --classes 100 --permutations 100000 --min-count 3 --min-coverage 6 --max-coverage 500 --chromosomes "2L 2R 3L 3R 4"
 # Input is a synchronized file with at least 4 populations; the first is the base-pop the other three the derived pop's
 # to speed it up use the script filter-synchronized before
  
=cut
