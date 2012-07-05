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
pod2usage(-verbose=>2) if $help;
pod2usage() unless -e $input;
pod2usage() unless $output;
pod2usage() unless -e $tehierfile;
pod2usage() unless $tehiertargetlevel;
pod2usage() unless $teinsertsitefile;


my $teh_resolver=get_te_hierarchy_resolver($tehierfile,$tehiertargetlevel);
my $te_inserts=load_te_inserts($teinsertsitefile);

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
print "Updating the polymorphism information\n";

my $poly_teinserts=[];
foreach my $tei (@$te_inserts)
{
    # frpes, fabs, rpres, rabs
    die "not existing" unless $tei; 
    my $polyi=$tei->subclone_polymorphism($tei->{tFpres},$tei->{tFabs},$tei->{tRpres},$tei->{tRabs});
    push @$poly_teinserts,$polyi;
}
print "Writing polymorphism data to file\n";
write_te_inserts($poly_teinserts,$output);
print "Potential TE-absence reads: $count_abs\n";
print "Potention TE-presence reads: $count_pres\n";

exit;



{
    package Utility;
    use strict;
    use warnings;
    
    sub get_presence_annotator
    {
        my $tei_h=shift;
        #$tei->{chr}{insdir}{pos} = []

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
            
            my $te_insertions=$tei_h->{$chr}{$insdir}{$pos};
            my $teids=$teresolver->($r->{chrmate});
            
            foreach my $tei (@$te_insertions)
            {
                if($tei->{teid} eq $teids)
                {
                    $tei->{"t".$insdir."pres"}++;
                }
            }
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
                my $te_insertions=$tei_h->{$chr1}{F}{$pos1};
                foreach my $te (@$te_insertions)
                {
                    $te->{tFabs}++;
                }
            }
            if(exists($tei_h->{$chr2}{R}{$pos2}))
            {
                my $te_insertions=$tei_h->{$chr2}{R}{$pos2};
                foreach my $te (@$te_insertions)
                {
                    $te->{tRabs}++;
                }
            }
            return;
        }
    }
    
    sub construct_insert_hash
    {
        my $inserts=shift;
        my $ih={};
        # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
        # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
        # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
        foreach my $tei (@$inserts)
        {
            my $chr=$tei->{chr};
            if($tei->{frstart})
            {
                $tei->{tFpres}=0; $tei->{tFabs}=0;
                my $start=$tei->{frstart};
                my $end=$tei->{frend};
                for my $i ($start..$end)
                {
                    $ih->{$chr}{F}{$i}||=[];
                    push @{$ih->{$chr}{F}{$i}},$tei;
                }
            }
            if($tei->{rrstart})
            {
                $tei->{tRpres}=0; $tei->{tRabs}=0;
                my $start=  $tei->{rrstart};
                my $end=    $tei->{rrend};
                for my $i ($start..$end)
                {
                    $ih->{$chr}{R}{$i}||=[];
                    push @{$ih->{$chr}{R}{$i}},$tei;
                }
            }
        }
        return $ih;
    }
    
}


    #"sam-file=s"            =>\$input,
    #"output=s"              =>\$output,
    #"min-map-qual=i"        =>\$minMapquality,
    #"te-hierarchy-file=s"   =>\$tehierfile,
    #"te-hierarchy-level=s"  =>\$tehiertargetlevel,
    #"te-insert-file=s"      =>\$teinsertsitefile,
    #"test"                  =>\$test,
    #"help"	            =>\$help

=head1 NAME

perl estimate-polymorphism.pl - Estimte the insertion frequencies for a given set of TE insertions

=head1 SYNOPSIS

 perl estimate-polymorphism.pl --sam-file pe_maped_reads.sam --te-insert-file te_insertions.txt --output te-insertion-polymorphism.txt --min-map-qual 15 --te-hierarchy-file te_hierarcy.txt --te-hierarchy-level family

=head1 OPTIONS

=over 4

=item B<--te-insert-file>

A file containing a list of TE insertions for which the insertion frequencies. Such a file is typically generated by <crosslink-te-sites.pl>; Mandatory

=item B<--sam-file>

A sam file containing the PE reads of a pooled population which have been mapped to the modified reference sequence

=item B<--output>

An output file containing the polymorphism information for every TE insertion of the input file

=item B<--min-map-qual>

the minimum mapping quality

=item B<--te-hierarchy-file>

a file containing a custom TE hierarchy

=item B<--te-hierarchy-level>

the level at the TE hierarchy at which the operations should be conducted

=item B<--help>

Display the help

=back

=head1 DETAILS


=cut
