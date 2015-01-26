use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use TEHierarchy;
use ParseSam;


my $input="";

my $te_ann="";
my $help=0;

GetOptions(
    "input=s"	        =>\$input,
    "te-annotation=s"   =>\$te_ann,
    "help"	        =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

my $an=get_TE_hierarchy($te_ann);


open my $ifh,"<",$input or die "Could not open sam file";


my($countsamefam,$countdiffam)=(0,0);

my $sp=get_basic_samparser();

while(my $l=<$ifh>)
{
    chomp $l;
    next if $l=~/^@/;
    my $s=$sp->($l);
    

    my $flag=$s->{flag};
    next if $flag & 0x0004;

    next unless $flag & 0x0001;
    next if $flag & 0x0008; # i am not interested in reads where the mate is unmapped

    
    
    my $cr=$s->{chr};
    my $cm=$s->{chrmate};
    if(exists($an->{$cr}) and exists($an->{$cm}))
    {
        my $rt=$an->{$cr}{family};
        my $mt=$an->{$cr}{family};
        if($rt eq $mt)
        {
            $countsamefam++;
        }
        else{
            $countdiffam++;
        }
    }
    
}

print $countsamefam;
print $countdiffam;
    


exit;