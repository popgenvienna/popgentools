#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my @inputFiles;
my $minMapQual=0;
my $help=0;

    
GetOptions(
           "input=s"        => \@inputFiles,

           "min_mapqual=i"  => \$minMapQual,
           "help"           => \$help
        );

pod2usage(-verbose=>2) if $help;
pod2usage(-verbose=>1,-message=>"You have to specify at least one input file") unless @inputFiles;


foreach (@inputFiles)
{
    pod2usage(-verbose=>1,-message=>"Input file $_ does not exist") unless -e $_;
}



#@SQ     SN:2L   LN:23011544
#read_no_mismatch        0       2L      229701  37      76M     *       0       0       GGTATTCGTATTAAAAAATATAAATTAATTATTGAATTGCATGTATAACTTGGCGGTTTTTGTATTATGAGGACGA    
#bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb    XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:76
 


my($countReads,$countMapped)=(0,0);
my $pairendflag=0;

my $mapQual=[];

# paired end statistics
my($mappedInProperPair, $fucked, $fuckedOtherContig,$fuckedSameStrand,$fuckedSamePosition)=(0,0,0,0,0,0,0);
my ($count_mpe_reads,$count_mse_reads,$count_spe_reads,$count_sse_reads)=(0,0,0,0);
my $fuckedMateUnmapped=0;
my $fuckedMatesOverlap=0;
my $fuckedImpropableDistance=0;
my $countUnmapped=0;


foreach my $file (@inputFiles)
{
    open my $ifh,"<",$file;
    while(my $e=<$ifh>)
    {
	chomp $e;
        next if $e=~/^@/; # discard the sam header
	my($readid,$flag,$refchr,$refpos,$mapqual,$cigar,$matechr,$matepos,$matedist,$read)=split /\s+/,$e;
	my $readlength=length($read);
        $countReads++;
	
	if($flag & 0x0001)
	{
	    $count_spe_reads++;
	}
	else
	{
	    $count_sse_reads++;
	}
        
	if($flag & 0x004)
	{
	    $countUnmapped++;
	    next;
	}
	else
	{
	    $countMapped++;    
	}
        
        
        #Calculate the mapping quality
        $mapQual->[$mapqual]||=0;
        $mapQual->[$mapqual]++;
        
        if($flag & 0x0001)
        {
	    $count_mpe_reads++;
            $pairendflag=1;
            
	    if($flag & 0x0002)
            {
                $mappedInProperPair++;
            }
            else 
            { # inproper pair
                $fucked++;
                
		if($flag & 0x0008) # only one read is mapped
		{
		    $fuckedMateUnmapped++;
		}
		else # both reads are mapped, but in an inproper pair
		{
		    if($matechr eq "=") # same contig
		    {
			if($matepos==$refpos)
			{
			    $fuckedSamePosition++;
			}
			else
			{	# not the same positon
			    
			    my $stranda=$flag & 0x0010?"R":"F";
			    my $strandb=$flag & 0x0020?"R":"F";

			    if($stranda eq $strandb)
			    {
				$fuckedSameStrand++;
			    }
			    else
			    {
				my $absdist=abs($matedist);
				if ($absdist<$readlength)
				{
				    $fuckedMatesOverlap++;
				}
				else
				{
				    $fuckedImpropableDistance++;
				}
			    }
			    
			}
			
		    }
		    else  # different contigs
		    {
			$fuckedOtherContig++;
		    }
		}

                
                
            }
            
            
        }
	else
	{
	    $count_mse_reads++;
	}
        
	
    }
    close $ifh;
}

print "Total number of reads: $countReads\n";
print "Single end reads: $count_sse_reads\n";
print "Paired end reads: $count_spe_reads\n\n";

print "Reads mapped: $countMapped\n";
print "Reads unmapped: $countUnmapped\n";
print "Mapped single end reads: $count_mse_reads\n";
print "Mapped paired end reads: $count_mpe_reads\n";

print "\n\nDistribution of mapping qualities (single+paired end)\n";
for(my $i=0; $i<@$mapQual; $i++)
{
    next unless $mapQual->[$i];
    print "$i\t$mapQual->[$i]\n";
}
print "\n\n";


print "Paired end statistic\n";
if($pairendflag)
{
    #($mappedInProperPair, $fucked, $fuckedOtherContig,$fuckedNegativeDist,$fuckedSameStrand,$fuckedSamePosition,$fuckedToFarAway)
    print "Mapped paired end reads: $count_mpe_reads\n";
    print "Mapped in proper pair: $mappedInProperPair\n";
    print "Mapped in inproper pair: $fucked\n";
    print "Mate is unmapped: $fuckedMateUnmapped\n";
    print "Mate is on different contig: $fuckedOtherContig\n";
    print "Mate is on the same strand: $fuckedSameStrand\n";
    print "Mates are identical (ie map to exactly the same position): $fuckedSamePosition\n";
    print "Mates overlap: $fuckedMatesOverlap\n";
    print "Impropable distance between mates: $fuckedImpropableDistance\n";
    
}


exit;

=head1 NAME

Sam-statistics.pl 

=head1 SYNOPSIS

 # Minimal argument call, specifying all required parameters.
 Sam-statistics.pl --input file1.sam --input file2.sam

=head1 OPTIONS

=over 4

=item B<--input>

The input file(s); has to be at least one sam file; Mandatory parameter.

=item B<--help>

Display the help pages

=back

=head1 DESCRIPTION

=head2 Input

A sam file as described here: http://samtools.sourceforge.net/


=head2 Output

Sam statistics


=head1 REQUIREMENTS

Perl 5.8 or higher

=head1 AUTHORS

Robert Kofler

=cut