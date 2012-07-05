#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use SAMPairReader;
my $input;
my $output;
my $outputduplicates;
my $help=0;
my $onlyproperpairs=0;
my $debug=0;    

GetOptions(
    
            "input=s"           => \$input,
            "output=s"          => \$output,
            "output-duplicates=s"=> \$outputduplicates,
            "proper-pairs"      => \$onlyproperpairs,
            "help"              => \$help,
            "debug"             =>\$debug
        ) or pod2usage(-verbose=>1,-message=>"invalid options");
pod2usage(-verbose=>2) if $help;
pod2usage(-verbose=>1,-message=>"You have to specify an input file") unless -e $input;

my $spr=SAMPairReader->new($input,100000);

open my $ifh,"<", $input or die "Could not open input file";
my $activechr="";
my $pastchromosomes={};
my $duplhash={};
my $stat={};
my $globalduplhash={};
my $pecount=0;
my $posduplicates=0;
while(my $pair=$spr->next())
{
    my($first,$second)=@$pair;

    
    if($onlyproperpairs)
    {
        next unless $first->{flag}  & 0x002;
        next unless $second->{flag} & 0x002;
    }
    next unless $first->{chr} eq $second->{chr};
    
    
    $activechr=$first->{chr} unless $activechr;
    # deal with chromosome changes; necessary to avoid excessive memory consumption
    if($activechr ne $first->{chr})
    {
        die "Samfile must be sorted by chromosome and position" if(exists($pastchromosomes->{$activechr}));
        $pastchromosomes->{$activechr}=1;
        Utility::updatestat($stat,$duplhash);
        Utility::updateglobalduplhahs($globalduplhash,$duplhash,$activechr) if $outputduplicates;
        $duplhash={};
        $activechr=$first->{chr};
    }
        
    #switch in case the first is the second
    ($first,$second)=($second,$first) if ($first->{start_s} > $second->{start_s});
    my $key = $first->{start_s}.":".$second->{end_s}; 
    $pecount++;
    $posduplicates++ if $first->{start_s} == $second->{start_s};
    $duplhash->{$key}++;
}
$spr->close;


# also update the last statistics
Utility::updatestat($stat,$duplhash);
Utility::updateglobalduplhahs($globalduplhash,$duplhash,$activechr) if $outputduplicates;
$duplhash={};

open my $ofh, ">",$output  or die "Could not open output file";
my $pairs = [ sort {$a->[0] <=> $b->[0] } map { [$_ , $stat->{$_}] } keys(%$stat) ];
my $posduplfraction=sprintf("%.4f",$posduplicates/$pecount);
print $ofh "Fragments considered\t$pecount\n";
print $ofh "Fragments where both mates have the same start postion\t$posduplicates\t$posduplfraction\n";

foreach my $p(@$pairs)
{
    my $dupllevel=$p->[0];
    my $count=$p->[1];
    my $fraction=$count*$dupllevel/$pecount;
    $fraction=sprintf("%.4f",$fraction);
    print $ofh "$dupllevel\t$count\t$fraction\n";
}
close $ofh;

# Finally parse the sam file again and print all sam entries having a duplicate signature (duplicate signatures are stored in globalduplhahsh)
# samfile, outputfile, onlyproperpairs, globalduplhash
Utility::print_globalduplicates($input,$outputduplicates,$onlyproperpairs,$globalduplhash) if $outputduplicates;

exit;

{
    package Utility;
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use SAMPairReader;

    sub updatestat
    {
        my $stat=shift;
        my $duplhash=shift;
        
        while(my($key,$count)=each(%$duplhash))
        {
            $stat->{$count}++;
        }
    }
    
    sub updateglobalduplhahs
    {
        my $globalduplhash=shift;
        my $duplhash=shift;
        my $activechr=shift;
        
        while(my($key,$count)=each(%$duplhash))
        {
            if($count>1)
            {
                $globalduplhash->{"$activechr:$key"}=1;
            }
        }
    }
    
    # samfile, outputfile, sp, onlyproperpairs, globalduplhash
    sub print_globalduplicates
    {
        my $samfile=shift;
        my $outputfile=shift;
        my $onlyproperpairs=shift;
        my $globalduplhash=shift;
            
        print "Print parsing samfile to extract the duplicates\n";
        open my $ofh, ">", $outputfile or die "Could not open output file";
        my $spr=SAMPairReader->new($samfile,100000);
        while(my $pair=$spr->next())
        {
            my($first,$second)=@$pair;
        
            
            if($onlyproperpairs)
            {
                next unless $first->{flag}  & 0x002;
                next unless $second->{flag} & 0x002;
            }
            next unless $first->{chr} eq $second->{chr};
            my $chr=$first->{chr};
            ($first,$second)=($second,$first) if ($first->{start_s} > $second->{start_s});
            my $key =$chr.":".$first->{start_s}.":".$second->{end_s}; 
            next unless(exists($globalduplhash->{$key}));
            
            print $ofh _formatsam($first)."\n";
            print $ofh _formatsam($second)."\n";
        }
    }
        
    sub _formatsam
    {

            my $sam=shift;
            my @entry=($sam->{readid},$sam->{flag},$sam->{chr},$sam->{start},$sam->{mq},$sam->{cigar},$sam->{chrmate},$sam->{posmate},$sam->{distance},$sam->{seq},$sam->{qual},$sam->{appendix});
            return join("\t",@entry); 
 
    }

    
}

=head1 NAME

PE-duplicate-statistics.pl 

=head1 SYNOPSIS

 # Minimal argument call
 PE-duplicate-statistics.pl --input input.sam

=head1 OPTIONS

=over 4

=item B<--input>

The input file(s); has to be SORTED (by chromosome and position) sam file. Mandatory parameter

=item B<--output>

Duplicate statistics; Mandatory parameter

=item B<--output-duplicates>

The duplicated paired end fragments in 'sam' format; Optional parameter

=item B<--proper-pairs>

Use only proper pairs for detecting duplicates

=item B<--help>

Display the help pages

=back

=head1 DESCRIPTION

=head2 Input

A sam file as described here: http://samtools.sourceforge.net/
Has to be SORTED (by chromosome and positions)

=head2 Output


=head1 REQUIREMENTS

Perl 5.8 or higher

=head1 AUTHORS

Robert Kofler

=cut