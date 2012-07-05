use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
our $verbose=1;


my $input="";

my $popoolationpath="";
my $help=0;
my $test=0;

my $minqual=20;
my $encoding="illumina";


#--gtf /Users/robertkofler/testfiles/3R-3entries.gtf --input /Users/robertkofler/testfiles/3R.pileup --output /Users/robertkofler/testfiles/output/3R-filtered.pileup

GetOptions(
    "input=s"           =>\$input,
    "popoolation-path=s"=>\$popoolationpath,
    "min-qual=i"        =>\$minqual,
    "fastq-type=s"      =>\$encoding,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";
pod2usage(-verbose=>2) if $help;
SubsTest::runTests() if $test;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $popoolationpath;


my $poplib=File::Spec->catfile($popoolationpath,"Modules");
eval("use lib '$poplib';");
eval("use Pileup;");


my $pp=get_basic_parser($encoding,$minqual);

open my $ifh,"<",$input or die "Could not open input file";

my $stat={};
while(my $l=<$ifh>)
{
    chomp $l;
    my $p=$pp->($l);
    my $cov=$p->{totcov};
    my $chr=$p->{chr};
    my $ac=[{a=>'A',c=>$p->{A}},
            {a=>'T',c=>$p->{T}},
            {a=>'C',c=>$p->{C}},
            {a=>'G',c=>$p->{G}},
            {a=>'N',c=>$p->{N}},
            {a=>'*',c=>$p->{del}}];
    $ac=[sort {$a->{c}<=>$b->{c}} @$ac ];
    
    my $mincount=0;
    # remove major allele
    pop @$ac;
    while(@$ac)
    {
        my $alc=shift @$ac;
        if($alc->{a} ne '*' and $alc->{a} ne 'N' and $alc->{c}>0)
        {
            $mincount=$alc->{c};
            last;
        }
    }
    $stat->{$cov}{$chr}[$mincount]++;
}


while(my($cov,$tmp)=each(%$stat))
{
    while(my($chr,$list)=each(%$tmp))
    {
        for(my $i=0; $i<@$list; $i++)
        {
            my $count=$list->[$i];
            next unless $count;
            print "$cov\t$chr\t$i\t$count\n";
        }
    }
}