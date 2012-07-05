use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use List::Util qw/min max/;
our $verbose=1;

my $input="";

my $popoolation2path="";
my $help=0;
my $test=0;




GetOptions(
    "input=s"           =>\$input,
    "popoolation2-path=s"=>\$popoolation2path,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";
pod2usage(-verbose=>2) if $help;
SubsTest::runTests() if $test;
pod2usage(-msg=>"Could not find sync file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Could not find path to PoPoolation2 file",-verbose=>1) unless -e $popoolation2path;
my $poplib=File::Spec->catfile($popoolation2path,"Modules");
eval("use lib '$poplib';");
eval("use Synchronized;");

my $sp=get_basic_syncparser();


open my $ifh, $input or die "Could not open input file";

my $costat=[];
my $mincov=[];
while (my $l=<$ifh>)
{
    chomp $l;
    my $p=$sp->($l);
    my $samp=$p->{samples};
    
    my @totcovs=();
    for(my $i=0; $i<@$samp; $i++)
    {
        my $cs=$samp->[$i];
        my $tc=$cs->{totcov};
        $costat->[$i][$tc]++;
        push @totcovs,$tc;
    }
    
    my $min=min(@totcovs);
    $mincov->[$min]++;
}

for(my $i=0; $i<@$costat; $i++)
{
    my $popnum=$i+1;
    my $tcov=$costat->[$i];
    
    my $sum=totcov($tcov);
    my $tcounter=0;
    for(my $k=0; $k<@$tcov; $k++)
    {
        my $cc=$tcov->[$k] || 0 ;
        $tcounter+=$cc;
        my $tfraction=sprintf("%.2f",$tcounter/$sum);
        print "$i\t$k\t$cc\t$tfraction\n";
    }
}

my $minsum=totcov($mincov);
my $mincounter=0;
for(my $i=0; $i<@$mincov; $i++)
{
    my $cc=$mincov->[$i] ||0;
    $mincounter+=$cc;
    my $tfraction=sprintf("%.2f",$mincounter/$minsum);
    print "min\t$i\t$cc\t$tfraction\n";
}



exit;


sub totcov
{
   my $covar=shift;
   
   
   my $sum=0;
   for(my $i=0; $i<@$covar; $i++)
   {
        my $count=0;
        $count=$covar->[$i] if $covar->[$i];
        $sum+=$count;
   }
    return $sum;
}