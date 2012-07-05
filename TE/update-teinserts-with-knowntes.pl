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
use TEInsertUtility;
use TEInsert;
use Utility;

my $teinsertfile;
my $tehierfile;
my $knownfile;
my $output;
my $tehiertargetlevel="family";
my $maxdist=300;
my $singlesiteshift=100;
my $help=0;
my $test=0;

GetOptions(
    "te-insertions=s"	    =>\$teinsertfile,
    "known=s"               =>\$knownfile,
    "output=s"              =>\$output,
    "te-hierarchy-file=s"   =>\$tehierfile,
    "te-hierarchy-level=s"  =>\$tehiertargetlevel,
    "max-distance=i"        =>\$maxdist,
    "single-site-shift=i"   =>\$singlesiteshift,
    "test"                  =>\$test,
    "help"	            =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

die "" unless $tehierfile;
die "" unless $tehiertargetlevel;
die "" unless $knownfile;
die "No TE insertion file provided" unless $teinsertfile;



my $teh_resolver=get_te_hierarchy_resolver($tehierfile,$tehiertargetlevel);
my $known=Utility::load_known_te_inserts($knownfile,$teh_resolver);
my $teinsertions=load_te_inserts($teinsertfile);


my $updated=[];
my $awaitfurtherprocessing=[];
my $statistics={};

my $knownidh={};


foreach my $ins (@$teinsertions)
{
 
    # INTERFACE:
    # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
    # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
    # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
    $statistics->{processed}++;
    my $chr=$ins->{chr};
    my $teid=$ins->{teid};
    my $ss=$ins->{sitesupport};

    

    my ($fwdknown,$revknown)=(undef,undef);
    if($ins->{frstart})
    {
        my $end=$ins->{frend};
        POS: for my $i($end..($end+$maxdist))
        {
            my $key="$chr:F:$teid:$i";
            if(exists($known->{$key}))
            {
                $fwdknown=$known->{$key};
                delete($known->{$key});
                last POS;
            }
        }
    }
    if($ins->{rrstart})
    {
        my $start = $ins->{rrstart};
        POS: for my $i (reverse(($start-$maxdist)..$start))
        {
            my $key="$chr:R:$teid:$i";
            if(exists($known->{$key}))
            {
                $revknown=$known->{$key};
                delete($known->{$key});
                last POS
            }
        }
    }
    
    if($ss eq "FR")
    {
        # TWO SITE support
        if($fwdknown and $revknown)
        {
            if($fwdknown eq $revknown)
            {
                $statistics->{correctpair}++;
                $ins->{fbid}=$fwdknown;
                push @$updated, $ins;
            }
            elsif($fwdknown ne $revknown)
            {
                # splitinsertions, mergedinsertions
                $statistics->{splitinsertions}++;
                my($fwd,$rev)=Utility::split_TEInsertion($ins,$singlesiteshift);
                $rev->{fbid}=$revknown;
                $fwd->{fbid}=$fwdknown;
                push @$awaitfurtherprocessing,$fwd;
                push @$awaitfurtherprocessing,$rev;
            }
            else
            {
                 die "impossible";
            }
        }
        elsif($fwdknown or $revknown)
        {
            my $fbid=$fwdknown?$fwdknown:$revknown;
            $ins->{fbid}=$fbid;
            push @$updated,$ins
        }
        else
        {
            push @$updated,$ins;
        }
    }
    elsif($ss eq "F" or $ss eq "R")
    {
        # ONE SITE supportingle sites
        if($revknown and $fwdknown)
        {
            die 'impossible';
        }
        elsif($revknown or $fwdknown)
        {
            my $fbid=$fwdknown?$fwdknown:$revknown;
            $ins->{fbid}=$fbid;
            push @$awaitfurtherprocessing,$ins;
        }
        else
        {
            push @$updated,$ins;
        }
    }
    else
    {
        die "unknown site support";
    }
    

}


# processing the singles which should be connected
my $unitihash={};
foreach my $tei (@$awaitfurtherprocessing)
{
    my $fbid=$tei->{fbid};
    die "No flybase id $fbid" unless $fbid;
    $unitihash->{$fbid}||=[];
    push @{$unitihash->{$fbid}},$tei;
}


while(my($fbid,$tes)=each(%$unitihash))
{
    if(@$tes==1)
    {
        push @$updated,$tes->[0];
    }
    elsif(@$tes==2)
    {
        my($fwd,$rev)=(undef,undef);
        $fwd=$tes->[0] if($tes->[0]{sitesupport} eq "F");
        $fwd=$tes->[1] if($tes->[1]{sitesupport} eq "F");
        $rev=$tes->[0] if($tes->[0]{sitesupport} eq "R");
        $rev=$tes->[1] if($tes->[1]{sitesupport} eq "R");
        unless($fwd and $rev)
        {
            die "Do not have both insertions for forward and reverse insertion";    
        }
        $statistics->{mergedinsertions}++;
        
        my $merged=Utility::merge_forwardandrev($fwd,$rev);
        $merged->{fbid}=$fbid;
        push @$updated,$merged;
    }
    else
    {
        die "to many elements for $fbid";
    }
}

# processed, bothidentified, correctpaired, wrongpaired, oneknown, zeroknown
# finalcount, knownidentified,knownwithone,knownwithtwo
$statistics->{finalcount}=@$updated;
my $fbidhash={};
foreach my $te (@$updated)
{
    my $fbid=$te->{fbid};
    next unless $fbid;
    my $sup=length($te->{sitesupport});
    $fbidhash->{$fbid}=$sup;
}
while(my($fbid,$count)=each(%$fbidhash))
{
    if($count>0)
    {
        $statistics->{knownidentified}++;
        $statistics->{knownwithone}++ if $count==1;
        $statistics->{knownwithtwo}++ if $count==2;
    }
}

write_te_inserts($updated,$output);

$statistics->{splitinsertions}=0 unless(exists($statistics->{splitinsertions}));
$statistics->{mergedinsertions}=0 unless(exists($statistics->{mergedinsertions}));
# splitinsertions, mergedinsertions
print "Processed count: $statistics->{processed}\n";
print "Final count: $statistics->{finalcount}\n";
print "Endstatistic\n";
print "Knownidentified: $statistics->{knownidentified}\n";
print "Known with one: $statistics->{knownwithone}\n";
print "Known with two: $statistics->{knownwithtwo}\n";

print "Merge split statistic\n";
print "Correct pair: $statistics->{correctpair}\n";
print "Splits: $statistics->{splitinsertions}\n";
print "Merges: $statistics->{mergedinsertions}\n";


exit;


{
    package Utility;
    use strict;
    use warnings;
    
    sub merge_forwardandrev
    {
        my $fwd     =shift;
        my $rev     =shift;
        die "can not merge sites of different chromosomes" if $fwd->{chr} ne $rev->{chr};
        my $inspos  =($fwd->{frend}+$rev->{rrstart})/2;
        my $teid    =$fwd->{teid};
        my $order   =$fwd->{order};
        my $fbid    =$fwd->{fbid};

        
        # INT: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
        my $merged=TEInsert->new($fwd->{chr},$inspos,$teid,$order,$fbid,"merge",
                                        $fwd->{frstart},$fwd->{frend},$fwd->{fpres},$fwd->{fabs},$fwd->{foverlap},
                                        $rev->{rrstart},$rev->{rrend},$rev->{rpres},$rev->{rabs},$rev->{roverlap});
        return $merged;
    }
    
    sub split_TEInsertion
    {
        my $ins=shift;
        my $singlesiteshift=shift;
        # INT: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
        my $fwd=TEInsert->new($ins->{chr},$ins->{frend}+$singlesiteshift,$ins->{teid},$ins->{order},$ins->{fbid},"split",$ins->{frstart},$ins->{frend},$ins->{fpres},$ins->{fabs},
                                $ins->{foverlap},undef,undef,undef,undef,undef);
        my $rev= TEInsert->new($ins->{chr},$ins->{rrstart}-$singlesiteshift,$ins->{teid},$ins->{order},$ins->{fbid},"split",undef,undef,undef,undef,
                                undef,$ins->{rrstart},$ins->{rrend},$ins->{rpres},$ins->{rabs},$ins->{roverlap});
        return ($fwd,$rev);
    }

    
    
    
    sub load_known_te_inserts
    {
        my $file=shift;
        my $resolver=shift;
        open my $ifh, "<", $file or die "Could not open input file";
        my $knownh={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            #2R	R	6867626	FBti0018865
            #2R	F	10972539	FBti0018866
            #2R	R	10979530	FBti0018866
            my($chr,$insdir,$pos,$fbid)=split /\t/,$l;
            my $teid=$resolver->($fbid);
            my $key="$chr:$insdir:$teid:$pos";
            $knownh->{$key}=$fbid;
        }
        return $knownh; #"$chr:$insdir:$teid:$pos"
    }
}

#GetOptions(
#    "te-insertions=s"	    =>\$teinsertfile,
#    "known=s"               =>\$knownfile,
#    "output=s"              =>\$output,
#    "te-hierarchy-file=s"   =>\$tehierfile,
#    "te-hierarchy-level=s"  =>\$tehiertargetlevel,
#    "max-distance=i"        =>\$maxdist,
#    "single-site-shift=i"   =>\$singlesiteshift,
#    "test"                  =>\$test,
#    "help"	            =>\$help


=head1 NAME

perl samro.pl - Merges two sam-files and crosslinks the paired end information (e.g.: when using bwa sw)

=head1 SYNOPSIS

 perl samro.pl --fq1 reads1.fastq --fq2 reads2.fastq --sam1 maped-reads1.sam --sam2 maped-reads2.sam --output merged.sam

=head1 OPTIONS

=over 4

=item B<--fq1>

The first read of a pair in fastq. Mandatory

=item B<--fq2>

The second read of a pair in fastq. Mandatory

=item B<--sam1>

the mapping results of the first read in sam. Mandatory

=item B<--sam2>

the mapping results of the second read in sam. Mandatory

=item B<--output>

The output file. [sam]  Mandatory.

=item B<--help>

Display the help

=back

=head1 DETAILS

=cut

