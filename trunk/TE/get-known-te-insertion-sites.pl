use strict;
use warnings;


my $file =shift;

open my $ifh, "<", $file or die "Could not open input file";


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
    print "$chr\tF\t$start\t$fbti\n";
    print "$chr\tR\t$end\t$fbti\n";


}