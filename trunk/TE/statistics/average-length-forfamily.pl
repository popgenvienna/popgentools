use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use TEHierarchy;
use Getopt::Long;
use Pod::Usage;


my $file =shift;
my $tehierfile=shift;

my $tr=get_te_hierarchy_translator($tehierfile,"insert","family");

open my $ifh, "<", $file or die "Could not open input file";

my $famlist={};


while(my $l=<$ifh>)
{
    next unless $l=~/^>/;
    chomp $l;
    #>FBti0018869 type=transposable_element; loc=2R:14463358..14469487; name=3S18{}853; dbxref=FlyBase_Annotation_IDs:TE18869,FlyBase:FBti0018869; MD5=99f36f1e2101d7e643ef5e072f2d231b; length=6130; release=r5.31; species=Dmel; 
    #>FBti0018871 type=transposable_element; loc=2R:complement(19801877..19809443); name=412{}880; dbxref=FlyBase_Annotation_IDs:TE18871,FlyBase:FBti0018871; MD5=09a65127014d3c89b32d233ebcf021d3; length=7567; release=r5.31; species=Dmel; 
    my($fbti)=$l=~/^>(\S+)/;
    my($loc)=$l=~/loc=([^;]+)/;
    my($chr,$start,$end)=$loc=~/^([^:]+):(?:complement\()?(\d+)\.\.(\d+)\)?$/;
    ($start,$end)=($end,$start) if $start>$end;
    my $leng=$end-$start+1;
    my $fam=$tr->{$fbti};
    $famlist->{$fam}||=[];
    push @{$famlist->{$fam}},$leng;
}


while(my($fam,$famcol)=each(%$famlist))
{
    $famcol=[sort{$a<=>$b} @$famcol];
    my $average=sprintf("%.1f",average_leng($famcol));
    my $min=$famcol->[0];
    my $max=$famcol->[-1];
    my $median=median_leng($famcol);
    print "$fam\t$min:$median:$average:$max\n";
}


exit;

sub median_leng
{
    my $list=shift;
    my $leng =@$list;
    if($leng%2)
    {
        my $idxlow=int($leng/2)-1;
        my $lower=$list->[$idxlow];
        my $upper=$list->[$idxlow+1];
        return ($lower+$upper)/2;
        
    }
    else
    {
        my $idx=int($leng/2);
        return $list->[$idx];
    }
}

sub average_leng
{
    my $list=shift;
    
    my $sum=0;
    $sum+=$_ foreach @$list;
    $sum/=scalar(@$list);
    return $sum
}