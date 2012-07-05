#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $inputfile;
my $outputFile;
my $feature="exon";
my $source="blat";
my $help=0;
my $ignoreMultiBest=1;


GetOptions(
           "input=s"        => \$inputfile,
           "output=s"       => \$outputFile,
           "feature=s"      => \$feature,
           "source=s"       => \$source,
           "print-multiple-hits" => sub{$ignoreMultiBest=0},
           "help"           => \$help
        ) or pod2usage(-msg=>"Invalid arguments",-verbose=>1);


pod2usage(-verbose=>2) if $help;

pod2usage(-msg=>"Could not find psl file",-verbose=>1) unless -e $inputfile;
pod2usage(-msg=>"Output file must be specified",-verbose=>1) unless $outputFile;

# psl; zero-based
# gtf: one-based, inclusive; end position is inclusve eg start_codon length 3: 380 ... 382 



open my $ofh, ">",$outputFile or die "Could not open output file";

###########################################################################################
########                     PSL SPECIFICATIONS                               #############
###########################################################################################
#
#
#psLayout version 3
#0          1       2   3       4           5       6       7       8       9           10      11      12      13                14    15      16      17      18              19          20
#match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
#        match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#2019    3       0       0       0       0       13      23954   +       gi|155369759|ref|NM_001101028.1|        2022    0       2022    chr12   41636750        2408426 2434402 14      129,160,161,108,120,116,170,163,191,180,106,144,207,67, 0,129,289,450,558,678,794,964,1127,1318,1498,1604,1748,1955,    2408426,2408841,2421866,2424191,2426442,2426708,2427245,2428829,2429838,2430599,2430904,2431141,2431434,2434335,
#363     0       0       0       0       0       2       2601    -       gi|147905018|ref|NM_001097514.1|        2222    0       363     chr13   135349118       96522098        96525062        3       143,88,132,     1859,2002,2090, 96522098,96522842,96524930,
#2201    8       0       0       1       1       13      15503   +       gi|178056503|ref|NM_001123076.1|        2244    34      2244    chr4    129318515       104497762       104515474       14      181,167,164,220,130,195,132,164,127,164,121,174,111,159,        34,215,383,547,767,897,1092,1224,1388,1515,1679,1800,1974,2085, 104497762,104498889,104502993,104506772,104507463,104508797,104509285,104510098,104511917,104512222,104513069,104514232,104514576,104515315,
#1244    6       0       0       0       0       4       4346    +       gi|165973425|ref|NM_001113706.1|        1256    0       1250    chr7    134056967       28517207        28522803        5       153,246,282,166,403,    0,153,399,681,847,      28517207,28520017,28520706,28521288,28522400,
#411     39      0       0       1       78      1       1265    +       gi|165973417|ref|NM_001113695.1|        1223    197     725     chr7    134056967       28691723        28693438        2       165,285,        197,440,        28691723,28693153,
#1174    23      0       0       2       26      5       11415   -       gi|165973417|ref|NM_001113695.1|        1223    0       1223    chr7    134056967       28594279        28606891        7       264,99,24,111,282,245,172,      0,265,364,388,499,781,1051,     28594279,28594543,28594926,28595430,28596548,28599216,28606719,
#906     32      0       0       4       113     5       4150    -       gi|165973417|ref|NM_001113695.1|        1223    172     1223    chr7    134056967       28575106        28580194        8       115,110,99,24,111,284,47,148,   0,154,265,364,388,499,855,903,  28575106,28575255,28575365,28575748,28576252,28577390,28579999,28580046,
#467     45      0       0       2       41      2       2177    -       gi|165973417|ref|NM_001113695.1|        1223    172     725     chr7    134056967       28561799        28564488        3       285,88,139,     498,788,912,    28561799,28564239,28564349,
#834     0       0       0       0       0       2       7050    +       gi|157073932|ref|NM_001103213.1|        850     1       835     chr4    129318515       128895665       128903549       3       72,289,473,     1,73,362,       128895665,128900554,128903076,
#834     0       0       0       0       0       2       7051    +       gi|157073932|ref|NM_001103213.1|        850     1       835     chr4    129318515       128804588       128812473       3       72,289,473,     1,73,362,       128804588,128809478,128812000,
#594     0       0       0       0       0       1       718     -       gi|148539981|ref|NM_214402.2|   609     15      609     chr14   148877259       48431109        48432421        2       411,183,        0,411,  48431109,48432238,
#37      0       0       0       1       7       1       1       +       gi|145279656|ref|NM_214014.2|   2335    472     516     chr8    92346407        77527344        77527382        2       20,17,  472,499,        77527344,77527365,
#47      0       0       0       2       7       2       2       +       gi|145279656|ref|NM_214014.2|   2335    474     528     chr5    83878508        15283745        15283794        4       5,19,15,8,      474,479,501,520,        15283745,15283751,15283770,15283786,
#40      2       0       0       0       0       0       0       +       gi|145279656|ref|NM_214014.2|   2335    474     516     chr12   41636750        1938031 1938073 1       42,     474,    1938031,
#32      2       0       0       0       0       0       0       +       gi|145279656|ref|NM_214014.2|   2335    469     503     chr12   41636750        10134498        10134532        1       34,     469,    10134498,
#53      1       0       0       2       6       2       24      +       gi|145279656|ref|NM_214014.2|   2335    456     516     chr11   65622336        11144771        11144849        4       5,7,31,11,      456,462,474,505,        11144771,11144776,11144789,11144838,
#33      2       0       0       0       0       0       0       +       gi|145279656|ref|NM_214014.2|   2335    478     513     chr1    264548797       255242227       255242262       1       35,     478,    255242227,
#37      2       0       0       0       0       0       0       -       gi|145279656|ref|NM_214014.2|   2335    476     515     chr9    109117252       104761684       104761723       1       39,     1820,   104761684,
#41      1       0       0       0       0       0       0       -       gi|145279656|ref|NM_214014.2|   2335    474     516     chr12   41636750        22509285        22509327        1       42,     1819,   22509285,
#40      2       0       0       0       0       0       0       -       gi|145279656|ref|NM_214014.2|   2335    474     516     chr12   41636750        2113368 2113410 1       42,     1819,   2113368,
#2158    53      0       0       1       109     3       95247   -       gi|145279656|ref|NM_214014.2|   2335    14      2334    chr11   65622336        13473522        13570980        4       121,1357,691,42,        1,122,1479,2279,        13473522,13474105,13570135,13570938,
#36      2       0       0       0       0       0       0       -       gi|145279656|ref|NM_214014.2|   2335    478     516     chr11   65622336        65538129        65538167        1       38,     1819,   65538129,
#153     6       0       0       4       24      4       31      +       gi|47523107|ref|NM_213875.1|    1033    658     841     chrX    90416161        43835460        43835650        6       16,4,26,72,7,34,        658,682,686,713,786,807,        43835460,43835481,43835486,43835512,43835591,43835616,
#
#############################################################################################################


    # 153     6       0       0       4       24      4       31      +       gi|47523107|ref|NM_213875.1|    1033    658     841     chrX    90416161        43835460        43835650        6       16,4,26,72,7,34,        658,682,686,713,786,807,        43835460,43835481,43835486,43835512,43835591,43835616,
    open my $ifh,"<",$inputfile or die "Could not open input file";
    my $start=0;
    my $best;
    HIT: while(<$ifh>)
    {
        # Check if the reads are already starting:
        if($_=~m/^\-+$/)
        {
            $start=1;
            next;
        }
        next unless $start;
        
        chomp;
        my @temp=split /\t/, $_;
        my $tLengths=$temp[18];
        my $tStart=$temp[20];
        my @Lengths=split /,/,$tLengths;
        my @Starts=split/,/,$tStart;
        
        my $hit = {
            match   =>  $temp[0],
            chr     =>  $temp[13],
            strand  =>  $temp[8],
            name    =>  $temp[9],
            lengths =>  \@Lengths,
            starts  =>  \@Starts
        };
        
        unless($best)
        {
            $best=[$hit];
            next HIT;
        }
        
        if($best->[0]{name} eq $hit->{name}) #only take the best hit for each
        {
            # query-id of the best is identical with current query-id
            if($hit->{match}==$best->[0]{match})
            {
                # the new match has the same score as the old match
                push @$best,$hit;
               
            }
            elsif($hit->{match} > $best->[0]{match})
            {
                 $best=[$hit];
            }
        }
        else
        {
            processBestHits($ofh,$best,$feature,$source,$ignoreMultiBest);
            $best=[$hit];
        }
    }
    
    processBestHits($ofh,$best,$feature,$source,$ignoreMultiBest);
    close $ifh;
    
exit(0);


sub processBestHits
{
    my $ofh=shift;
    my $besthits=shift;
    my $feature=shift;
    my $source=shift;
    my $ignoreMultiHits=shift;
    
    my $hitcount=@$besthits;
    
    if($hitcount==1)
    {
        _parseHit($ofh,$_,$feature,$source) foreach(@$besthits);
    }
    elsif($hitcount> 1 and $ignoreMultiHits)
    {
        # don't do anything
    }
    elsif($hitcount>1 and not $ignoreMultiHits)
    {
        _parseHit($ofh,$_,$feature,$source) foreach(@$besthits);
    }
    else
    {
        die "impossible combination";
    }
}

sub _parseHit
{
    my $ofh     = shift;
    my $hit     = shift;
    my $feature = shift;
    my $source  = shift;
    
    my @starts=@{$hit->{starts}};
    my @lengths=@{$hit->{lengths}};
    my $chr=$hit->{chr};
    my $strand=$hit->{strand};
    my $name=$hit->{name};

    die "length vector has to have the same size as the start vector" unless scalar(@starts)==scalar(@lengths);

    for(my $i=0;$i<@starts;$i++)
    {
        my $start=$starts[$i];
        my $leng=$lengths[$i];
        
        _gtfHit($ofh,$chr,$feature,$source,$start,$leng,$strand,$name);
    }
}


sub _gtfHit
{


    my $ofh=shift;
    my $chrom=shift;
    my $feature=shift;
    my $source=shift;
    my $chromStart=shift;
    my $length=shift;
    my $strand=shift;
    my $name=shift;
    
    # gtf is 1-based whereas psl is zero-based
    $chromStart++;
    
    # gtf is inclusive
    my $chromEnd=$chromStart+$length-1;
    
    my $score=".";
    print $ofh "$chrom\t$source\t$feature\t$chromStart\t$chromEnd\t$score\t$strand\t.\tgene_id \"$name\"; transcript_id \"$name\";\n";
    
}

           #"input=s"        => \$inputfile,
           #"output=s"       => \$outputFile,
           #"feature=s"      => \$feature,
           #"source=s"       => \$source,
           #"print-multiple-hits" => sub{$ignoreMultiBest=0},
           #"help"           => \$help
           
=head1 NAME

perl psl2gtf.pl - A script which converts selects the BEST blat hit and converts it into a gtf file

=head1 SYNOPSIS

perl psl2gtf.pl --input input.psl --output best-blat-hit.gtf

=head1 OPTIONS

=over 4

=item B<--input>

The input file;  a psl-file;  Mandatory.

=item B<--output>

A gtf file as specified here: http://mblab.wustl.edu/GTF2.html

=item B<--feature>

The feature which is processed; default=exon

=item B<--source>

The source of the gtf file; default=blat

=item B<--print-multiple-hits>

flag; enable printing of multiple best hits; per default this option is disabled, i.e: multiple best hits will be ignored

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

a psl file

=head2 Output

a gtf file
  
=cut
