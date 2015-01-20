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
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

my $an=get_TE_hierarchy($te_ann);
my $teh_resolver=get_te_hierarchy_resolver($te_ann,"family");
open my $ifh,"<",$input or die "Could not open sam file";


my $families=["1360","17.6","1731","297","3S18","412","accord","accord2","aurora-element","baggins","Bari1","Bari2","blood","BS","BS3","BS4","Burdock","Circe","copia","Cr1a","diver","diver2","Dm88","Doc","Doc2-element","Doc3-element","Doc4-element","F-element","FB","flea","frogger","Fw2","Fw3","G-element","G2","G3","G4","G5","G5A","G6","G7","GATE","gtwin","gypsy","gypsy10","gypsy11","gypsy12","gypsy2","gypsy3","gypsy4","gypsy5","gypsy6","gypsy7","gypsy8","gypsy9","HB","Helena","HeT-A","HMS-Beagle","HMS-Beagle2","hobo","hopper","hopper2","I-element","Idefix","INE-1","invader1","invader2","invader3","invader4","invader5","invader6","Ivk","jockey","jockey2","Juan","looper1","Mariner","mariner2","Max-element","McClintock","mdg1","mdg3","micropia","NOF","opus","Osvaldo","P-element","pogo","Porto1","Q-element","Quasimodo","R1-2","R1A1-element","R2-element","roo","rooA","rover","Rt1a","Rt1b","Rt1c","S-element","S2","springer","Stalker","Stalker2","Stalker3","Stalker4","Tabor","TAHRE","Tc1","Tc1-2","Tc3","Tirant","Tom1","transib1","transib2","transib3","transib4","Transpac","X-element","ZAM"];
my $euchr={"X"=>1, "2L"=>1, "2R"=>1, "3L"=>1, "3R"=>1, "4"=>1};





## read stats
my $mapedpres=0;
my $mapedboth=0;

my $sp=get_basic_samparser();
my $allchrstat={};
my $euchrstat={};

while(my $l=<$ifh>)
{
    chomp $l;
    next if $l=~/^@/;
    my $s=$sp->($l);
    my $flag=$s->{flag};
    next if $flag & 0x0004;
    next unless $flag & 0x0001;
    next if $flag & 0x0008; # i am not interested in reads where the mate is unmapped
    $mapedboth++;

    my $cr=$s->{chr};
    my $cm=$s->{chrmate};
    
    if( not exists($an->{$cr}) and exists($an->{$cm}))
       {
        $mapedpres++;
        my $fam=$teh_resolver->($cm);
        # cr maps to chromosome and mate to a TE
        $allchrstat->{$fam}||=0;
        $allchrstat->{$fam}++;
        if(exists($euchr->{$cr}))
        {
            $euchrstat->{$fam}||=0;
            $euchrstat->{$fam}++;
        }
        
       }
}
    
foreach my $fam (@$families)
{
    my $allc=0;
    my $euc=0;
    $allc=$allchrstat->{$fam} if(exists($allchrstat->{$fam}));
    $euc=$euchrstat->{$fam} if(exists($euchrstat->{$fam}));
    print "$fam\t$allc\t$euc\n";
}
print "maped\t$mapedboth\n";
print "mapedpres\t$mapedpres\n";



        
exit;