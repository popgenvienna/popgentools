#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use List::Util qw[min max];
use TEHierarchy;
use ParseSam;
use TEInsertUtility;
use TESamReader;

my $makenoise=100000;
my $input;
my $tehierfile;
my $tehiertargetlevel;
my $teinsertsitefile;
my $readlength=75;
my $output="";
my $minMapquality=20;
my $help=0;
my $test=0;

GetOptions(
    "sam-file=s"            =>\$input,
    "output=s"              =>\$output,
    "min-map-qual=i"        =>\$minMapquality,
    "te-hierarchy-file=s"   =>\$tehierfile,
    "te-hierarchy-level=s"  =>\$tehiertargetlevel,
    "te-insert-file=s"      =>\$teinsertsitefile,
    "test"                  =>\$test,
    "help"	            =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);
pod2usage() unless -e $input;
pod2usage() unless $output;
pod2usage() unless -e $tehierfile;
pod2usage() unless $tehiertargetlevel;
pod2usage() unless $teinsertsitefile;


my $teh_resolver=get_te_hierarchy_resolver($tehierfile,$tehiertargetlevel);
my $te_inserts=load_te_inserts($teinsertsitefile);

# remove overlapping TE inserts
$te_inserts=get_nonoverlapping_te_inserts($te_inserts);
# create an TE insertion hash $tei->{chr}{insdir}{pos} = $teinsert
my $tei_hash=Utility::construct_insert_hash($te_inserts);

my $samparser=get_te_samparser($teh_resolver);
my $sr=TESamReader->new($input,$samparser,$makenoise);

my($abs_annotator,$pre_annotator)=(Utility::get_absence_annotator($tei_hash,$minMapquality,$teh_resolver),Utility::get_presence_annotator($tei_hash,$minMapquality,$teh_resolver));
print "Start reading the Sam file\n";
my($count_abs,$count_pres)=(0,0);
while(my $s=$sr->next())
{
    if($s->{mode} eq "abs")
    {
        $count_abs++;
        $abs_annotator->($s);
    }
    elsif(  $s->{mode} eq "pre")
    {
        $count_pres++;
        $pre_annotator->($s);
    }
    else
    {
        die "Mode not supported";
    }
}
print "Finished reading the sam file\n";
print "Potential TE-absence reads: $count_abs\n";
print "Potention TE-presence reads: $count_pres\n";
print "\nStart calculating the Transposon pi for each element\n";
open my $ofh, ">",$output or die "Could not open output file";
foreach my $tei (@$te_inserts)
{
            # entry: chr, insdir, teid, count, start, end, count_pre, count_abs;
    my $abs=$tei->{count_abs};
    my $pre=$tei->{count_pre};
    my $sum=$abs+$pre;
    my $te_pi="na";
    $te_pi=1-($abs/$sum)**2 - ($pre/$sum)**2 if $sum;
    print $ofh "$tei->{chr}\t$tei->{insdir}\t$tei->{teid}\t$tei->{start}\t$tei->{end}\t$te_pi\t$tei->{count_pre}\t$tei->{count_abs}\n";
    
}



exit;



{
    package Utility;
    use strict;
    use warnings;
    
    sub get_presence_annotator
    {
        my $tei_h=shift;
        #$tei->{chr}{insdir}{pos} 
        # entry: chr, insdir, teid, count, start, end, count_pre, count_abs
        my $minqual=shift;
        my $teresolver=shift;
        return sub
        {
            my $ins=shift;
            my $r=$ins->{r1};
            # readid, flag, chr, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
            # start, start_s, end, end_s, te_ins_read, te_ins_direction, pp
            my $chr=$r->{chr};
            my $insdir=$r->{te_ins_direction};
            die "Insertion direction not specified" if ($insdir eq "-");
            die "Error: $chr is not a reference contig; It has been identified as TE" if($teresolver->($chr));
            die "Error: $r->{chrmate} is not a TE" unless($teresolver->($r->{chrmate}));
            my $pos=$insdir eq "F"? $r->{start_s}:$r->{end_s};
            
            #discard if mapping quality is not sufficient;
            return if $r->{mq}<$minqual;
            
            # discard if not in the range hash
            return unless exists($tei_h->{$chr}{$insdir}{$pos});
            
            my $te=$tei_h->{$chr}{$insdir}{$pos};
            my $teids=$teresolver->($r->{chrmate});
            
            return if($teids ne $te->{teid});
            
            # increase the presence count;
            $te->{count_pre}++;
            return;
        }
    }
    
    sub get_absence_annotator
    {
        my $tei_h=shift;
        #$tei->{chr}{insdir}{pos} 
        # entry: chr, insdir, teid, count, start, end, count_pre, count_abs
        my $minqual=shift;
        my $teresolver=shift;
        return sub
        {
            my $ins=shift;
            my $r1=$ins->{r1};
            my $r2=$ins->{r2};
            my $chr1=$r1->{chr};
            my $chr2=$r2->{chr};
            die "Chromosomes of pair are not identical $chr1 and $chr2\n"unless $chr1 eq $chr2;
            die "Read 1 maps to a transposon" if ($teresolver->($chr1));
            die "Read 2 maps to a transposon" if ($teresolver->($chr2));
            
            # check if the two reads are not overlapping
            return unless($r1->{end_s}< $r2->{start_s});
            
            
            return if $r1->{mq}< $minqual;
            return if $r2->{mq}< $minqual;
            # readid, flag, chr, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
            # start, start_s, end, end_s, te_ins_read, te_ins_direction, pp
            die "Strand of the insert is wrong; has to be fwd" if($r1->{flag} & 0x0010);
            die "Strand of the insert is wrong; has to be ref" unless($r2->{flag} &0x0010);
            

            my $pos1=$r1->{start_s};
            my $pos2=$r2->{end_s};
            
            
            if(exists($tei_h->{$chr1}{F}{$pos1}))
            {
                my $te=$tei_h->{$chr1}{F}{$pos1};
                $te->{count_abs}++;
            }
            if(exists($tei_h->{$chr2}{R}{$pos2}))
            {
                my $te=$tei_h->{$chr2}{R}{$pos2};
                $te->{count_abs}++;
            }
            return;
        }
    }
    
    sub construct_insert_hash
    {
        my $inserts=shift;
        my $ih={};
        # chr, insdir, teid, count, start, end, count_pre, count_abs
        foreach my $tei (@$inserts)
        {
            my $chr=$tei->{chr};
            my $insdir=$tei->{insdir};
            my $start=$tei->{start};
            my $end=$tei->{end};
            for my $i ($start..$end)
            {
                $ih->{$chr}{$insdir}{$i}=$tei;
            }
        }
        return $ih;
    }
    
}

