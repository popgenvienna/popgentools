use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use SamSlider;


my $input="";
my $output="";
my $windowsize  = 10_000;
my $stepsize    = 10_000;
my $measure="mean";
my $help=0;
my $test=0;

GetOptions(
    "input=s"               =>\$input,
    "output=s"              =>\$output,
    "window-size=i"         =>\$windowsize,
    "step-size=i"           =>\$stepsize,
    "median"                => sub{$measure="median";},
    "help"                  =>\$help,
    "test"                  =>\$test
    ) or die "Could not parse arguments";

pod2usage(-verbose=>2) if $help;
PairTest::runTests() if $test;

pod2usage(-msg=>"Could not find input file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;

open my $ofh, ">$output" or die "Could not open output file";

my $reader = SamSlider->new($input,$windowsize,$stepsize);

my $method;
if ($measure eq "mean")
{
    $method=\&calculate_mean;
}
elsif($measure eq "median")
{
    $method=\&calculate_median;
}
else
{
    die "invalid measure";
}

while(my $window=$reader->nextWindow())
{
        my $chr=$window->{chr};
        my $pos=$window->{middle};
        my $win=$window->{window};
        my $data=$window->{data};
        
        next unless @$data;
        # check if all reads are paired end reads in sequence
        foreach my $r(@$data)
        {
            warn "Sam-file contains reads which are not paired in sequencing\n" unless $r->{p}; # paired in sequence
        }
         
        my $dist = get_distances($data);
        my $dist_count=@$dist;
        
        my $meas=$method->($dist);
        $meas=sprintf("%.2f",$meas) if $meas ne "na";
        
        print $ofh "$chr\t$pos\t$dist_count\t-\t$meas\n";
}
exit;


sub calculate_median
{
    my $distances=shift;
    my @dist=sort {$a<=>$b} @$distances;
    
    my $index = int(scalar(@dist)/2);
    my $val=$dist[$index];
    return "na" unless $val;
    return $val;
}


sub calculate_mean
{
    my $distances=shift;
    my $sum=0;
    my $count=0;
    foreach(@$distances)
    {
        $count++;
        $sum+=$_;
    }
    return "na" unless $count;
    return $sum/$count;
}

sub get_distances
{
    my $data=shift;
    my $distances=[];
    
    foreach my $sam (@$data)
    {
        next unless $sam->{P}; # mapped in a proper pair
        next unless $sam->{posmate}>$sam->{pos};
        
        my $cigar=$sam->{cigar};
        my $leng=alignmentlengthFromCigar($cigar);
        
        my $dist=$sam->{posmate}-($sam->{pos}+$leng-1);

        push @$distances,$dist;
    }
    return $distances;    
}


    sub alignmentlengthFromCigar
    {
        my $cigar=shift;
        my (@e)=$cigar=~/(\d+[MSDIN])/g;
        
        
        my $leng=0;
        
        for(my $i=0; $i<@e; $i++)
        {
            my $en=$e[$i];
            
            if($en=~/^(\d+)M/)
            {
                $leng+=$1;
            }
            elsif($en=~/^(\d+)S/)
            {
                
            }
            elsif($en=~/^(\d+)I/)
            {
 
            }
            elsif($en=~/^(\d+)[ND]/)
            {
                $leng+=$1;
            }
            else
            {
                die "not allowed";
            }
        }
        return $leng;
    }


{
    package PairTest;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use SamSlider;
    use Test;
    use Test_SamSlider;
    
    sub runTests
    {
        test_SamSlider();
        test_mean();
        test_median();
        exit;
    }
    
    sub test_mean
    {
        my $mean;
        $mean=main::calculate_mean([1,2]);
        is($mean,"1.5","calculate mean; mean is ok");
        $mean=main::calculate_mean([3]);
        is($mean,"3","calculate mean; mean is ok");
        $mean=main::calculate_mean([3,4,5]);
        is($mean,"4","calculate mean; mean is ok");        
    }
    
    sub test_median
    {
        my $med;
        $med=main::calculate_median([1]);
        is($med,"1","calculate median; median is ok");
        $med=main::calculate_median([1,2]);
        is($med,"2","calculate median; median is ok");
        $med=main::calculate_median([1,2,3]);
        is($med,"2","calculate median; median is ok");
        $med=main::calculate_median([1,2,3,4]);
        is($med,"3","calculate median; median is ok");        
        
    }
}



=head1 NAME

perl pair-distance.pl - A script which calculates the average distance of paired end reads along contigs (chromosomes) using a sliding window approach.

=head1 SYNOPSIS

perl pair-distance.pl --input input.sam --output output.txt --window-size 10000 --step-size 10000

=head1 OPTIONS

=over 4

=item B<--input>

The input file; A SORTED! sam file; Mandatory

=item B<--output>

The output file.  Mandatory.

=item B<--median>

Flag, calculate median; Per default the script calculates the mean distance in a given window. Using this flag the median will be calculated instead.

=item B<--window-size>

The size of the sliding window. default=10000

=item B<--step-size>

the size of one sliding window step. If this number is equal to the C<--window-size> the sliding window will be non overlapping. default=10000

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

Usually paired end reads mapped to contigs in form of a  SORTED! sam file.
For a detailed description see: http://samtools.sourceforge.net/SAM1.pdf

Sorting can be done using the samtools!

=head2 Output

For every sliding window the mean or median distance between the pairs will be displayed

 2L	1730000	20	-	230.10
 2L	1740000	53	-	229.89
 2L	1750000	41	-	231.10
 2L	1760000	39	-	232.10
 2L	1770000	39	-	242.20
 
 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: number of proper pairs in the window, i.e.: number of pairs used for estimating col5
 col 4: empty (-); 
 col 5: inner distance between pairs (mean, or median)

=head2 Attention

When using sampe make sure to increase the allowed distance between paired end reads otherwise an artificial limit will be impossed and the distribution may erronously appear quite even.

This script calculates the inner distance between pairs;
It uses the start position of the first read to determine the position; Soft shaded ends (S in cigar) are not considered. Insertions (I) Deletions (D) and padding (N) are considered.

=cut


