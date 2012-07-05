use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use Pileup;

my $noise=10000;

my $inputfile=shift;
my $mincount=shift;
die "Provide a minimum count" unless $mincount>0;

my $pp=get_basic_mpileupparser("illumina",20);


open my $ifh,"<", $inputfile or die "could not open input file";
my $stat=[];
my $linecount=1;
while(my $line=<$ifh>)
{
    chomp $line;
    my $pe=$pp->($line);
    my $refc=$pe->{refc};
    my $entries=$pe->{entries};
    for(my $i=0; $i<@$entries; $i++)
    {
        my $e=$entries->[$i];
        my $ar=[{a=>"A",c=>$e->{A}},{a=>"T",c=>$e->{T}},{a=>"C",c=>$e->{C}},{a=>"G",c=>$e->{G}}];
        $ar=[sort {$a->{c}<=>$b->{c}} @$ar]; # sort ascending
        pop @$ar; # remov major allele
        
        my $countminor=0;
        my $allele="";
        foreach my $s (@$ar)
        {
            if($s->{c}>0)
            {
                $countminor=$s->{c};
                $allele=$s->{a};
                last;
            }
        }
        next unless $countminor;
        next if $refc eq $allele;
        next if $refc eq "N";
        next if $allele eq "N";
        
        if($countminor <= $mincount)
        {
            my $key=$refc.$allele;
            $stat->[$i]{$key}+=$countminor;
        }
    }
    warn "Processed $linecount pileupentries" unless($linecount % $noise);
    $linecount++;
    
}


for (my $i=0; $i<=@$stat; $i++)
{
    my $stath=$stat->[$i];
    next unless $stath;
    while(my($key,$count)=each(%$stath))
    {
        my $num=$i+1;
        print "$num\t$key\t$count\n";
    }
}

