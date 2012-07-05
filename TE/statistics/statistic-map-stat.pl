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
my $reads=0;
my $maped=0;
my $mapped_in_pair=0;
my ($pp)=0;
my $weird_pp=0;
my ($pp_chr,$pp_te)=(0,0);
my ($map_dif_cont)=0;
my ($map_dif_te,$map_dif_chr,$map_reads_te_insertion)=(0,0,0);

## pair stats
my ($pair_fwd_insertion,$pair_rev_insertion)=(0,0);
my ($pair_useful,$pair_mq20)=(0,0);

my $chimchr={};
my $tecounter={};
my $chrcounter={};
my $chrmate_mq=[];
my $temate_mq=[];

my $tecountcontrol={};

my $sp=get_basic_samparser();

while(my $l=<$ifh>)
{
    chomp $l;
    next if $l=~/^@/;
    my $s=$sp->($l);
    
    $reads++;
    my $flag=$s->{flag};
    next if $flag & 0x0004;
    $maped++;
    next unless $flag & 0x0001;
    next if $flag & 0x0008; # i am not interested in reads where the mate is unmapped
    $mapped_in_pair++;
    
    
    my $cr=$s->{chr};
    my $cm=$s->{chrmate};
    
    
    if($cm eq "=")
    {
        if($flag & 0x0002)
        {
            # proper pair and on the same contig......
            warn "$l\n" unless $cm eq "=";
        
            $pp++;
            if(exists($an->{$cr}) )
            {
                $pp_te++;
            }
            else
            {
                $pp_chr++;
            }
        }
    }
    elsif($cr ne $cm )
    {
        die "must not happen $cr $cm" if $cm eq "*" or $cm eq "=";
        $weird_pp++ if($flag & 0x0002);
        
        $map_dif_cont++;
        if(exists($an->{$cr}) and exists($an->{$cm}))
        {
            $map_dif_te++;
        }
        elsif(not exists($an->{$cr}) and not exists($an->{$cm}))
        {
            $map_dif_chr++;
            my @ar=($cr,$cm);
            @ar =sort {$a cmp $b} @ar;
            my $key="$ar[0]:$ar[1]";
            $chimchr->{$key}++;
            
            
        }
        elsif( (not exists($an->{$cr}) and exists($an->{$cm})) or (exists($an->{$cr}) and not exists($an->{$cm})) )
        {
            $map_reads_te_insertion++;
            my $mq=$s->{mq};
            
            if(not exists($an->{$cr}))
            {
                # $cr should be chromosome;
                $pair_useful++;
                $tecounter->{$cm}||=0;
                $chrcounter->{$cr}||=0;
                $tecounter->{$cm}++;
                $chrcounter->{$cr}++;
                
                $chrmate_mq->[$mq]++;
                
                if($flag & 0x0010)
                {
                    $pair_rev_insertion++;
                }
                else
                {
                    $pair_fwd_insertion++;
                }
                
            }
            else
            {
                $temate_mq->[$mq]++;
            }
        }
        else
        {
            $tecountcontrol->{$cr}||=0;
            $tecountcontrol->{$cr}++;
            die "strange: $l\n";
        }
        
    }
}
    
    print "stat\t[reads] all\t$reads\n";
    print "stat\t[reads] mapped\t$maped\n";
    print "stat\t[reads] both mates mapped\t$mapped_in_pair\n";
    print "stat\t[reads] mapped in proper pair (pp)\t$pp\n";
    print "stat\t[reads] weird proper pair (different contigs)\t$weird_pp\n";
    print "stat\t[reads] pp in chr\t$pp_chr\n";
    print "stat\t[reads] pp in te\t$pp_te\n";
    print "stat\t[reads] mates are on two different contigs\t$map_dif_cont\n";
    print "stat\t[reads] mates are on two different chromosomes\t$map_dif_chr\n";
    print "stat\t[reads] mates are on two different TEs\t$map_dif_te\n";
    print "stat\t[reads] one mate is on chromosome the other on TE\t$map_reads_te_insertion\n";
    print "stat\t[fragments] chimeric between chr and TE\t$pair_useful\n";
    print "stat\t[fragments] chimeric (TE - chr) and chr-read has mq>=20\t$pair_mq20\n";
    print "stat\t[fragments] forward insertions\t$pair_fwd_insertion\n";
    print "stat\t[fragments] reverse insertions\t$pair_rev_insertion\n";
    
for my $i (0..scalar(@$chrmate_mq)-1)
{
    next unless $chrmate_mq->[$i];
    print "chrmate_mq\t[fragments]\t$i\t$chrmate_mq->[$i]\n";
}
for my $i (0..scalar(@$temate_mq)-1)
{
    next unless $temate_mq->[$i];
    print "temate_mq\t[fragments]\t$i\t$temate_mq->[$i]\n";
}
while(my($te,$count)=each(%$tecounter))
{
    print "te_ins\t[fragments]\t$te\t$count\n";
}
while(my($te,$count)=each(%$tecountcontrol))
{
    print "te_ins_cont\t[fragments]\t$te\t$count\n";
}
while(my($k,$count)=each(%$chimchr))
{
    my ($k1,$k2)=split /:/,$k;
    print "chim_chr\t[reads]\t$k1\t$k2\t$count\n";
}


exit;