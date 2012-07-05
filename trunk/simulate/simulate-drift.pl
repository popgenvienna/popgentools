use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use List::Util qw[min max];
use Getopt::Long;
use Pod::Usage;

my $generations =0;
my $popSize     =0;
my $simulations =0;
my $avcov       =0;
my $help        =0;



GetOptions(
    "generations=s"	        =>\$generations,
    "pop-size=i"                =>\$popSize,
    "simulations=i"             =>\$simulations,
    "average-coverage=i"        =>\$avcov,
    "help"                      =>\$help            
    );

pod2usage(-msg=>"provide --generations") unless $generations;
pod2usage(-msg=>"provide --pop-size") unless $popSize;
pod2usage(-msg=>"provide --simulations") unless $simulations;
pod2usage(-msg=>"provide --average-coverage") unless $avcov;
pod2usage(-verbose=>2) if $help;

my $allelefreqdistribution=Utility::_get_allele_freq_distribution();
my $coveragedistribution=[];

my ($gen1,$gen2)=split /:/,$generations;


my $fsts=[];
for my $i (1..$simulations)
{
    my $targetfreq=$allelefreqdistribution->[int(rand(scalar(@$allelefreqdistribution)))];
    my $basepop=Utility::get_basepopulation($targetfreq,$popSize);
    my $dpop1=Utility::drift_population($basepop,$gen1);
    my $dpop2=Utility::drift_population($basepop,$gen2);
    
    
    my($ab,$ad1,$ad2)=(Utility::measure_allelefreq($basepop),Utility::measure_allelefreq($dpop1),Utility::measure_allelefreq($dpop2));
    my($fst1,$fst2,$fst12)=(Utility::calculateFst($ab,$ad1,$popSize),Utility::calculateFst($ab,$ad2,$popSize),Utility::calculateFst($ad1,$ad2,$popSize));
    
    print "$i\tall\tbase:1\t$fst1\t$targetfreq\t$ab->{A}\t$ab->{C}\t$ad1->{A}\t$ad1->{C}\n";
    print "$i\tall\tbase:2\t$fst2\t$targetfreq\t$ab->{A}\t$ab->{C}\t$ad2->{A}\t$ad2->{C}\n";
    print "$i\tall\t1:2\t$fst12\t$targetfreq\t$ad1->{A}\t$ad1->{C}\t$ad2->{A}\t$ad2->{C}\n";
    
    my($bpcov,$dpopcov1,$dpopcov2)=(Utility::random_draw($basepop,$avcov),Utility::random_draw($dpop1,$avcov),Utility::random_draw($dpop2,$avcov));
    my($abcov,$adcov1,$adcov2)=(Utility::measure_allelefreq($bpcov),Utility::measure_allelefreq($dpopcov1),Utility::measure_allelefreq($dpopcov2));
    my($fstcov1,$fstcov2,$fstcov12)=(Utility::calculateFst($abcov,$adcov1,$popSize),Utility::calculateFst($abcov,$adcov2,$popSize),Utility::calculateFst($adcov1,$adcov2,$popSize));
    
    print "$i\trd\tbase:1\t$fstcov1\t$targetfreq\t$abcov->{A}\t$abcov->{C}\t$adcov1->{A}\t$adcov1->{C}\n";
    print "$i\trd\tbase:2\t$fstcov2\t$targetfreq\t$abcov->{A}\t$abcov->{C}\t$adcov2->{A}\t$adcov2->{C}\n";
    print "$i\trd\t1:2\t$fstcov12\t$targetfreq\t$adcov1->{A}\t$adcov1->{C}\t$adcov2->{A}\t$adcov2->{C}\n";

# i 1:2 Fst A1 C1 A2 C2    
}









exit;

{
    package Utility;
    use strict;
    use warnings;
    
    # MINOR ALLELE = A
    # MAJOR ALLELE = C
    
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
    
    sub calculateFst
    {
        # calculte the Fst; provided are the allele counts for two populations
        my $pop1=shift;
        my $pop2=shift;
        my $popSize=shift;
        
        my $totpop=
        {
        A=>($pop1->{A}+$pop2->{A})/2,
        C=>($pop1->{C}+$pop2->{C})/2,
        cov=>($pop1->{cov}+$pop2->{cov})/2
        };
        
        my $pip1=_pi($pop1,$popSize);
        my $pip2=_pi($pop2,$popSize);
        my $pitot=_pi($totpop,$popSize);
        
        my $av=($pip1+$pip2)/2;
        
        my $fst=0;
        $fst=(($pitot-$av)/$pitot) if $pitot;
        return $fst;
    }
    
    sub _pi
    {
        my $snp = shift;
        my $popSize=shift;
        my $M=$snp->{cov};
        my $pi_snp=1;
        $pi_snp-=(($snp->{A})/$M)**2;
        $pi_snp-=(($snp->{C})/$M)**2;
    
        #this term is indepent of the correction
        $pi_snp*=($M/($M-1));
        $pi_snp*=($popSize/($popSize-1));
        return $pi_snp;
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



    sub _get_allele_freq_distribution
    {
        my $alleledistri=[[0.025,119],[0.075,59],[0.125,32],[0.175,21],[0.225,17],[0.275,16],[0.325,14],[0.375,13],[0.425,12],[0.475,12]];
        my $afd=[];
        foreach my $al (@$alleledistri)
        {
            my $af=$al->[0];
            my $count=$al->[1];
            for (1..$count)
            {
                push @$afd,$af;
            }
            
        }
        return $afd;
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
    
    
    #sub _randompairChromosomes
    #{
    #    my $pop=shift;
    #    my $pairs=[];
    #    while(@$pop)
    #    {
    #        my $aal=splice(@$pop,int(rand(scalar(@$pop))),1);
    #        my $bal=splice(@$pop,int(rand(scalar(@$pop))),1);
    #        push @$pairs,[$aal,$bal];
    #    }
    #    return $pairs;
    #}
    
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

=head1 NAME

simulate-drift.pl - simulates drift and calculates the Fst which may be expected due to genetic drift

=head1 SYNOPSIS
    
perl simulate-drift.pl --pop-size 500 --average-coverage 100 --generations 15:20 --simulations 1000

=cut

