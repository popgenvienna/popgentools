use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use VarianceExactCorrection;
use FormatSNP;
use Pileup;
our $verbose=1;

my $pileupfile="";
my $output="";
my $minCount=2;
my $fastqtype="illumina";
my $minQual=20;
my $poolSize=0;
my $help=0;
my $test=0;
my $minCoverage=4;
my $maxCoverage=400;
my $minCoveredFraction=0.6;
my $featureLength=0;
my $measure="";
my $snpfile;
my $snpsonly=0;


GetOptions(
    "measure=s"         =>\$measure,
    "input=s"           =>\$pileupfile,
    "output=s"          =>\$output,
    "snp-output=s"      =>\$snpfile,
    "fastq-type=s"      =>\$fastqtype,
    "min-count=i"       =>\$minCount,
    "min-qual=i"        =>\$minQual,
    "feature-length=i"  =>\$featureLength,
    "pool-size=i"       =>\$poolSize,
    "min-coverage=i"    =>\$minCoverage,
    "max-coverage=i"    =>\$maxCoverage,
    "min-covered-fraction=f"=>\$minCoveredFraction,
    "test"              =>\$test,
    "snpsonly"          =>\$snpsonly,
    "help"              =>\$help
) or die "Invalid arguments";


pod2usage(-verbose=>2) if $help;
VarTest::runTests() if $test;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $pileupfile;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;
pod2usage(-msg=>"You have to provide the length of the feature",-verbose=>1) unless $featureLength;
pod2usage(-msg=>"Pool size not provided",-verbose=>1) unless $poolSize;
pod2usage(-msg=>"Min count not provided",-verbose=>1) unless $minCount;
pod2usage(-msg=>"Min quality not valid. Has to be between 0 and 40",-verbose=>1) if $minQual<0 || $minQual > 40;
pod2usage(-msg=>"The minimum coverage hast to be at least two times the minimum count",-verbose=>1) unless $minCoverage >= (2*$minCount);
pod2usage(-msg=>"Measure not provided",-verbose=>1) unless $measure;
pod2usage(-msg=>"A SNP output file is mandatory when only a SNP output has been requested",-verbose=>) if($snpsonly and not $snpfile);



my $varianceCalculator=VarianceExactCorrection->new($poolSize,$minCount);

# qualencoding,mincount,mincov,maxcov,minqual
my $pp=get_pileup_parser($fastqtype,$minCount,$minCoverage,$maxCoverage,$minQual);

my $snpwriter=get_SNPwriter($snpfile) if $snpfile;


my ($snpcount,$covercount,$data)=Utility::read_file($pileupfile,$pp);

die "The feature length provided by the user can not be correct; Feature length must be larger or equal as the number of entries in the pileup file" if $covercount>$featureLength;


open my $ofh, ">",$output or die "Could not open output file $output\n";


    my $coveredFraction=$covercount/$featureLength;
    
    my $snps=[];
    foreach(@$data)
    {
        if ($_->{ispuresnp})
        {
            push @$snps,$_;
            $snpwriter->($_) if $snpwriter;
        }
    }

    unless($snpsonly)
    {
        my $meas=0;    
        $meas=$varianceCalculator->calculate_measure($measure,$snps,$covercount);
        $meas=sprintf("%.9f",$meas);
        
        $meas="na" if $coveredFraction < $minCoveredFraction;
        
        $coveredFraction=sprintf("%.3f",$coveredFraction);
        print $ofh "$snpcount\t$coveredFraction\t$meas\n";
    }


close $ofh;
exit;



{
    package Utility;
    use warnings;
    use strict;
    use FindBin qw($RealBin);
    use lib "$RealBin/../Modules";
    use VarianceExactCorrection;
    use Pileup;
    
    sub read_file
    {
        my $file=shift;
        my $pp=shift;
        
        open my $ifh, "<", $file or die "could not open file handle";
        
        my $data=[];
        while(my $l=<$ifh>)
        {
            chomp $l;
            my $p=$pp->($l);
            push @$data,$p;
            # $line, $mincount, $mincoverage,$minqual,$maxcoverage;
        }
        my $snps=0;
        my $covered=0;
        
        foreach my $p (@$data)
        {
            $snps++ if $p->{ispuresnp};
            $covered++ if $p->{iscov};
        }
        
        return ($snps,$covered,$data);
    }
    

}




=head1 NAME

perl Variance-in-file.pl - A script which calculates the requested population genetics measure (pi, theta, d) for a pileup file

=head1 SYNOPSIS

perl variance-slider.pl --measure pi --input input.pileup --output output.file --pool-size 500 --min-count 2 --min-coverage 4 ---max-coverage 400 --feature-length 650

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the pileup format. A pooled population sequenced and mapped to the reference genome. Finally the mapping results have been converted to sam output format.
Using the samtools the sam format can be easily converted into the pileup format.  Mandatory.

=item B<--output>

The output file;  Mandatory.

=item B<--snp-output>

If provided, this file will contain the polymorphic sites which have been used for this analysis; default=empty

=item B<--feature-length>

The length of the feature; Mandatory

=item B<--measure>

Currently, "pi", "theta" and "D" is supported. This stands for Tajima's Pi, Watterson's Theta, Tajima's D respectively; Mandatory

=item B<--pool-size>

The size of the pool which has been sequenced. e.g.: 500; Mandatory

=item B<--fastq-type>
The encoding of the quality characters; Must either be 'sanger' or 'illumina'; 

 Using the notation suggested by Cock et al (2009) the following applies:
 'sanger'   = fastq-sanger: phred encoding; offset of 33
 'solexa'   = fastq-solexa: -> NOT SUPPORTED
 'illumina' = fastq-illumina: phred encoding: offset of 64
 
 See also:
 Cock et al (2009) The Sanger FASTQ file format for sequecnes with quality socres,
 and the Solexa/Illumina FASTQ variants; 

default=illumina
 
=item B<--min-count>

The minimum count of the minor allele. This is important for the identification of SNPs; default=2

=item B<--min-coverage>

The minimum coverage of a site. Sites with a lower coverage will not be considered (for SNP identification and coverage estimation); default=4

=item B<--max-coverage>

the maximum coverage; used for SNP identification, the coverage in ALL populations has to be lower or equal to this threshold, otherwise no SNP will be called. default=1000000

=item B<--min-covered-fraction>

the minimum fraction of a window being between min-coverage and max-coverage in ALL populations; float; default=0.6

=item B<--min-qual>

The minimum quality; Alleles with a quality lower than this threshold will not be considered (for SNP identification and coverage estimation); default=20


=item B<--snpsonly>

Only output the SNPs; No Variance information will be calculated

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

A pileup file as described here: http://samtools.sourceforge.net/pileup.shtml; example:

 2L	90131	N	11	AaAAAaaAaAA	[aUQ_a`^_\Z
 2L	90132	N	11	AaAAAaaAaAA	_bYQ_^aaT^b
 2L	90133	N	11	A$aAAAaaAaAA	_b[Xaaa__Ua
 2L	90134	N	10	tTTTttTtTT	_`aaa_a[aa
 2L	90135	N	10	a$TAAaaAaAA	aZ^a`ba`\_
 2L	90136	N	9	TTTttTtTT	`aaaaaWaa
 2L	90137	N	9	GGGggGgGG	``aaaaQaa
 2L	90138	N	9	T$TTttTtTT	[U\`a\T^_
 2L	90139	N	8	TTttTtTT	``aaaU_a
 2L	90140	N	9	CCccCcCC^FC	[aaba`aaa

=head2 Output

The output will be as in the followin example:

 21	1.000	0.001365804


 col 1: SNPs found in the feature; the SNPs of all partial features - eg. exons of a gene - are summed up
 col 2: fraction of the window covered by a sufficient number of reads. Suficient means higher than min-coverage and lower than max-coverage; all partial features are considered
 col 3: population genetics estimator (pi, theta, d); this is the weighted mean of all partial features; weighted by length of the partial feature

=head2 SNP Output

 2R      53601   A       89      86      1       0       2       0
 2R      53642   T       105     0       81      24      0       0
 2R      53663   A       117     98      19      0       0       0

 The individual tab-delimited entries are in the following format:
 col 1: chromosome ID (contig)
 col 2: position in chromosome
 col 3: reference character
 col 4: coverage
 col 5: counts of A's
 col 6: counts of T's
 col 7: counts of C's
 col 8: counts of G's
 col 9: counts of N's


=head2 Technical Details

Note: This script has been specifically been developed for Ram's webservice; I can not think of any other useful applications; Don't use this script if you do not exactly know what to do
 
 
=cut
