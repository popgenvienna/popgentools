{
package SyncTEUtility;
use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use TEInsert;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(format_tesync load_tesync);



sub format_tesync
{
    my $s=shift;
    my $chr=$s->{chr};
    my $pos=$s->{pos};
    my $support=$s->{sitesupport};
    my $fam=$s->{fam};
    my $ord=$s->{ord};
    my $fbid=$s->{fbid};
    my $comment=$s->{comment};
    my $overlap=$s->{overlap};
    my $pop=$s->{pop};
    # chr, pos, sitesupport, fam, ord, fbid, comment, overlap, pop=>abs,pres
    $ord||="-";
    $fbid||="-";
    $comment||="-";
    # 2L      3004472 FR      INE-1   1       TIR     FBti0063010     ncorrdist=119   3004251 3004341 1       10      10      0       0       3004603 3004702 1       9       9       0       0
    my $str="$chr\t$pos\t$support\t$fam\t$ord\t$fbid\t$comment\t$overlap\t";
    
    my $tofo=[];
    foreach my $temp (@$pop)
    {
        my $tstr="$temp->{pres}:$temp->{abs}";
        push @$tofo,$tstr;
    }
    
    my $prob=join("\t",@$tofo);
    $str.=$prob;
    return $str;
}

sub load_tesync
{
    my $te_insert_file=shift;
    open my $ifh, "<", $te_insert_file or die "Could not open input file";
    my $syncte_inserts=[];
    
    while(my $l=<$ifh>)
    {
        chomp $l;
        next if $l=~/^#/;
        my $parsed=_parse_tesync($l);
        push @$syncte_inserts,$parsed;
    }
    return $syncte_inserts;
}

sub _parse_tesync
{
    my $line=shift;
    my @a=split /\t/,$line;
    my $e={};
    $e->{chr}=shift @a;
    $e->{pos}=shift @a;
    $e->{sitesupport}=shift @a;
    $e->{teid}=shift @a;
    $e->{ord}=shift @a;
    $e->{fbid}=shift @a;
    $e->{comment}=shift @a;
    $e->{overlap}=shift @a;
    
    my $pops=[];
    my $index=1;
    foreach my $pop (@a)
    {
        my ($pres,$abs)=split /:/,$pop;
        push @$pops,
        {
            pres=>$pres,
            abs=>$abs,
            cov=>$pres+$abs,
            index=>$index
        };
        $index++;
    }
    $e->{pops}=$pops;
    
    # chr, pos, sitesupport, teid, ord, fbid, comment, overlap, pops {pres, abs, cov, index}
    return $e;
}

}
1;
