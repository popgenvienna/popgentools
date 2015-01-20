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
my $families=["1360","17.6","1731","297","3S18","412","accord","accord2","aurora-element","baggins","Bari1","Bari2","blood","BS","BS3","BS4","Burdock","Circe","copia","Cr1a","diver","diver2","Dm88","Doc","Doc2-element","Doc3-element","Doc4-element","F-element","FB","flea","frogger","Fw2","Fw3","G-element","G2","G3","G4","G5","G5A","G6","G7","GATE","gtwin","gypsy","gypsy10","gypsy11","gypsy12","gypsy2","gypsy3","gypsy4","gypsy5","gypsy6","gypsy7","gypsy8","gypsy9","HB","Helena","HeT-A","HMS-Beagle","HMS-Beagle2","hobo","hopper","hopper2","I-element","Idefix","INE-1","invader1","invader2","invader3","invader4","invader5","invader6","Ivk","jockey","jockey2","Juan","looper1","Mariner","mariner2","Max-element","McClintock","mdg1","mdg3","micropia","NOF","opus","Osvaldo","P-element","pogo","Porto1","Q-element","Quasimodo","R1-2","R1A1-element","R2-element","roo","rooA","rover","Rt1a","Rt1b","Rt1c","S-element","S2","springer","Stalker","Stalker2","Stalker3","Stalker4","Tabor","TAHRE","Tc1","Tc1-2","Tc3","Tirant","Tom1","transib1","transib2","transib3","transib4","Transpac","X-element","ZAM"];
my $euchr=["X", "2L", "2R", "3L", "3R", "4"];


open my $ifh,"<",$input or die "Could not open sam file";



## read stats

my $maped=0;
my $mapedte=0;

my $sp=get_basic_samparser();
my $famstat={};

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
        my $fam=$an->{$chr}{family};
        $famstat->{$fam}++;
        $mapedte++;
    }
   
}

foreach  my $fam (@$families)
{
    my $count=0;
    $count=$famstat->{$fam} if(exists($famstat->{$fam}));
    my $fraction=0.0;
    $fraction=(100.0*$count)/$mapedte if $mapedte>0;
    print "$fam\t$count\t$fraction\n";
    
}

print "Mapped\t$maped\n";
print "Mapped TE\t$mapedte\n";
    

exit;