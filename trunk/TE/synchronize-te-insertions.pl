#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use Getopt::Long;
use Pod::Usage;
use Utility;
use TEInsertUtility;
use SyncTEUtility;


my @inputFiles;
my $output;

GetOptions(
    "input=s"  =>\@inputFiles,
    "output=s" =>\$output
) or die "Invalid parameters";

pod2usage() unless @inputFiles;
pod2usage() unless $output;


my $tereader= [ map {get_te_insertreader($_)} @inputFiles];


open my $ofh, ">", $output or die "Could not open output file handle\n";
while(1)
{
    my $nextreads= [map {$_->()} @$tereader];
    last unless $nextreads->[0];
    # chr, pos, sitesupport, fam, ord, comment, overlap, pop
    my $p = _create_sync_te($nextreads);
    my $formated=format_tesync($p);
    print $ofh $formated."\n";
    
}



exit;

sub _create_sync_te
{
    my $nextreads=shift;
    
    my $chr;
    my $pos;
    my $support;
    my $fam;
    my $ord;
    my $fbid;
    my $comment;
    my $overlap;
    my $presabs=[];
    
#2L      1527101 R       Juan    0       non-LTR -       -       -       -       -       -       -       -       -       1527201 1527296 0       27      0       27      0
#2L      1527618 F       jockey  0       non-LTR -       -       1527455 1527518 0       23      0       23      0       -       -       -       -       -       
#2L      3004472 FR      INE-1   1       TIR     FBti0063010     ncorrdist=119   3004251 3004341 1       10      10      0       0       3004603 3004702 1       9       9       0       0

    foreach my $r (@$nextreads)
    {
            # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
            # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
            # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
            $chr    =$r->{chr} unless $chr;
            $pos    =$r->{inspos} unless $pos;
            $support=$r->{sitesupport} unless $support;
            $fam    =$r->{teid} unless $fam;
            $chr eq $r->{chr} or die "Chromosome is not fitting $chr vs $r->{chr}\n";
            $pos == $r->{inspos} or die "Position is not fitting $pos vs $r->{inspos}\n";
            $support eq $r->{sitesupport} or die "Site support is not fitting $support vs $r->{sitesupport}\n";
            $fam eq $r->{teid} or die "Family is not fitting $fam vs $r->{teid}\n";
            
            $ord= $r->{order};
            $fbid=$r->{fbid};
            $comment=$r->{comment};
            $overlap=0;
            $overlap=1 if $r->{foverlap};
            $overlap=1 if $r->{roverlap};
            
            my $prescount=0;
            my $abscount=0;
            
            $prescount+=    $r->{fpres}     if $r->{fpres};
            $abscount+=     $r->{fabs}      if $r->{fabs};
            $prescount+=    $r->{rpres}     if $r->{rpres};
            $abscount+=     $r->{rabs}      if $r->{rabs};
            
            push @$presabs,{
                pres=>$prescount,
                abs=>$abscount
            };
    }
    return
    {
        chr=>$chr,
        pos=>$pos,
        sitesupport=>$support,
        fam=>$fam,
        ord=>$ord,
        fbid=>$fbid,
        comment=>$comment,
        overlap=>$overlap,
        pop=>$presabs
    };
    # chr, pos, sitesupport, fam, ord, fbid, comment, overlap, pop=>abs,pres
}

