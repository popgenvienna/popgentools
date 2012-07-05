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



## read stats

my $maped=0;
my $mapedte=0;

my $sp=get_basic_samparser();

while(my $l=<$ifh>)
{
    chomp $l;
    next if $l=~/^@/;
    my $s=$sp->($l);
    
    my $flag=$s->{flag};
    next if $flag & 0x0004;
    $maped++;

    my $chr=$s->{chr};
    if(exists($an->{$chr}))
    {
        $mapedte++;
    }
   
}

print "Mapped\t$maped\n";
print "Mapped TE\t$mapedte\n";
    

exit;