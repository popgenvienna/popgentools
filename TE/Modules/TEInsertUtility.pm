{
package TEInsertUtility;
use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use TEInsert;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(load_te_inserts write_te_inserts get_te_insertreader);

sub get_te_insertreader
{
    my $te_insert_file=shift;
    open my $ifh ,"<", $te_insert_file or die "Could not open input file";
    return sub
    {
        my $line=<$ifh>;
        return undef unless $line;
        chomp $line;
        return _parse_teinsert($line);
    }
}

sub load_te_inserts
{
    my $te_insert_file=shift;
    open my $ifh, "<", $te_insert_file or die "Could not open input file";
    my $te_inserts=[];
    while(my $l=<$ifh>)
    {
        chomp $l;
        my $tei=_parse_teinsert($l);
        push @$te_inserts,$tei;
    }
    return $te_inserts;
}

sub write_te_inserts
{
    my $teinserts=shift;
    my $file=shift;
    
    open my $ofh, ">",$file or die "Could not open output file";
    foreach my $tei (@$teinserts)
    {
        my $formated=$tei->tostring;
        print $ofh $formated."\n";
    }
    close $ofh;
}



sub _parse_teinsert
{
    my $l=shift;
        # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
        # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
        # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
    my($chr,$inspos,$sitesupport,$teid,$popfreq,$order,$fbid,$comment,$frstart,$frend,$fpopfreq,$fcov,$fpres,$fabs,$foverlap,$rrstart,$rrend,$rpopfreq,$rcov,$rpres,$rabs,$roverlap) = split /\t/,$l;
    $order eq "-"   and $order =undef;
    $fbid eq "-"    and $fbid=undef;
    $comment eq "-" and $comment=undef;
    
    $frstart eq "-" and $frstart=undef;
    $frend eq "-"   and $frend=undef;
    $fpres eq "-"   and $fpres=undef;
    $fabs eq "-"    and $fabs=undef;
    $foverlap eq "-"and $foverlap=undef;
    
    $rrstart eq "-" and $rrstart=undef,
    $rrend eq "-"   and $rrend=undef;
    $rpres eq "-"   and $rpres=undef;
    $rabs eq "-"    and $rabs=undef;
    $roverlap eq "-"and $roverlap=undef;
    
    die "Invalid TE Insertion $l\n" unless $chr;
    die "Invalid TE Insertion $l\n" unless $inspos;
    die "Invalid TE Insertion $l\n" unless $sitesupport;
    die "Invalid TE Insertion $l\n" unless $teid;
    die "Invalid TE Insertion $l\n" unless ($frstart or $rrstart);
    
    return TEInsert->new($chr,$inspos,$teid,$order,$fbid,$comment, $frstart,$frend,$fpres,$fabs,$foverlap,$rrstart,$rrend,$rpres,$rabs,$roverlap);
    # INT: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
    
    
}


}

1;
