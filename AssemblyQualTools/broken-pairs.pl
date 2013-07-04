use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use SamSlider;
use ParseSam;

my $input="";
my $output="";
my $windowsize  = 10_000;
my $stepsize    = 10_000;
my $help=0;
my $test=0;

GetOptions(
    "input=s"               =>\$input,
    "output=s"              =>\$output,
    "window-size=i"         =>\$windowsize,
    "step-size=i"           =>\$stepsize,
    "help"                  =>\$help,
    "test"                  =>\$test
    ) or die "Could not parse arguments";

pod2usage(-verbose=>2) if $help;
runTests() if $test;
pod2usage(-msg=>"Could not find input file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;

open my $ofh, ">$output" or die "Could not open output file";

my $reader = SamSlider->new($input,$windowsize,$stepsize);

while(my $window=$reader->nextWindow())
{
        my $chr=$window->{chr};
        my $pos=$window->{middle};
        my $win=$window->{window};
        my $data=$window->{data};
        
        next unless @$data;
            # ss..sam-statisitcs
            # mall=>$countall,   # mapped all
            # mpp=>$mappedInProperPair,
            # mse=>$countSingleEnd, # single ended reads mapped as single end
            # f_mipp=>$fucked, # fucked: mapped in improper pair
            # f_mum=>$fuckedMateUnmapped,
            # f_sp=>$fuckedSamePosition,
            # f_ss=>$fuckedSameStrand,
            # f_mo=>$fuckedMatesOverlap,
            # f_id=>$fuckedImpropableDistance,
            # f_oc=>$fuckedOtherContig  
        my $ss=calculate_sam_statistics($data);
        print $ofh "$chr\t$pos\t$ss->{mall}\t$ss->{mpp}\t$ss->{f_mipp}\t$ss->{f_mum}\t$ss->{f_oc}\t$ss->{f_ss}\t$ss->{f_id}\t$ss->{f_sp}\t$ss->{f_mo}\n";
}
exit;

=head1 NAME

perl broken-pairs.pl - A script which calculates the number of reads not mapped as proper pair using a sliding window approach

=head1 SYNOPSIS

perl broken-pairs.pl --input input.sam --output output.txt --window-size 10000 --step-size 10000

=head1 OPTIONS

=over 4

=item B<--input>

The input file; A SORTED! sam file; Mandatory

=item B<--output>

The output file. Mandatory.

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

For every sliding window a couple of parameters will be calculated
 
 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: number of reads mapped in the window (starting position)
 col 4: number of reads mapped in a proper pair 
 col 5: number of reads not mapped in a proper pair; not proper pair = npp
 col 6: npp; mate is unmapped
 col 7: npp; mate mapped to a different contig
 col 8: npp; mate mapped to the same strand
 col 9: npp; distance between the paired end reads is impropable
 col 10: npp; mate has exactly the same position (the paired end reads are identical)
 col 11: npp; the two paired end reads are overlapping (small fragment)

Note that col4 and col5 must not necessarily add up to col3, as col3 also contains reads which are not paired in sequence (single end reads)

=head2 Reformating for Visualise-output.pl using awk

Visualise output requires col1, col2 and col5 where col1 has to be the reference chromosome, col2 the position and col5 the measure which should be displayed.
Following some examples how the output of this script may be reformated using awk

 a.) Display the reads where the mate maps to a different contig
 awk 'BEGIN{OFS="\t"}{print $1,$2,"-","-",$7}' test.bp
 
 b.) Display the reads which are not maped in a proper pair normalised by all reads mapped to the window
 awk 'BEGIN{OFS="\t"}{print $1,$2,"-","-",$7/$3}' test.bp


=cut