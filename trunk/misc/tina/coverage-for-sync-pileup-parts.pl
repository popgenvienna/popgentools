#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use POSIX;

my $FILE_NAME;
my $PARTS_FILE_NAME;
my $INSERTIONS=2;
my @POP_COLUMNS;
my $HELP;

GetOptions(
	"input=s"=>\$FILE_NAME,
	"prettybase-summary=s"=>\$PARTS_FILE_NAME,
	"type=i"=>\$INSERTIONS, # 0 for SNPs, 1 for insertions, 2 for both types
	"population=s"=>\@POP_COLUMNS,
	'help'=>\$HELP,
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

if (!defined($FILE_NAME) or !defined($PARTS_FILE_NAME) or !defined(@POP_COLUMNS)){die "Invalid arguments, use 'perl $0 --help'."}

@POP_COLUMNS = (split(/,/, join(',', @POP_COLUMNS)));

my %genomicParts=();
open partsFileHandle, "<", $PARTS_FILE_NAME;
while (my $line = <partsFileHandle>){
	next if $line =~ m/^#/;
	chomp($line);
	my ($position, $type, $sig, $length) = split "\t", $line;
		
	$genomicParts{$position}{length}=$length;
	$genomicParts{$position}{type}=$type;
	$genomicParts{$position}{sig}=$sig;
}
close partsFileHandle;

#count total length of region of interest

my $totalLength=0;
foreach my $chr (keys %genomicParts){
  foreach my $pos (keys %{$genomicParts{$chr}}){
    if ($INSERTIONS == 0 and $genomicParts{$chr}{$pos}{type} eq 'SNP'){
	$totalLength+=1;
    }elsif($INSERTIONS == 1 and $genomicParts{$chr}{$pos}{type} eq 'INS'){
	$totalLength+=$genomicParts{$chr}{$pos}{length};
    }elsif($INSERTIONS == 2){
      $totalLength+=$genomicParts{$chr}{$pos}{length};
    }
  }
}

open fileHandle, "<", $FILE_NAME or die "Could not open file ".$FILE_NAME;
my $numberOfPopulations = scalar @POP_COLUMNS;
my @tmpSums = ();

#my @tmpCounts = ();

#print Dumper(\%genomicParts);
#print Dumper($totalLength);
my $in_insertion = 0;
my $lasts = 0;
my $ins = "";

while (my $line = <fileHandle>){
	next if $line =~ m/^#/;
	chomp($line);
	
	my ($chromosome, $position, $reference, @parts) = split "\t", $line;
	
	#print(exists($genomicParts{$chromosome}{$position}));
	#print($in_insertion);
      
      
	if (exists($genomicParts{$chromosome}{$position}) or ($in_insertion)){
		
		$ins = $genomicParts{$chromosome}{$position}{type};

		if (exists($genomicParts{$chromosome}{$position}) and ($genomicParts{$chromosome}{$position}{length} > 1)){
			$in_insertion=1;	
			$lasts = $genomicParts{$chromosome}{$position}{length};			
		}
		
		#print($INSERTIONS);
		#print("ins",$ins,"\n");
		
		if ((($INSERTIONS==1) and ($ins eq "INS")) or (($INSERTIONS==0) and exists($genomicParts{$chromosome}{$position}) and ($ins eq "SNP")) or ($INSERTIONS==2)){
			
			#print("in\n");
			foreach my $column (@POP_COLUMNS){
			
				my $actCov;
			
				if ($parts[$column-4] eq "-"){
					$actCov=0;	
				}else{
					my @pop = split ":", $parts[$column-4];
					$actCov = $pop[0]+$pop[1]+$pop[2]+$pop[3];
				}
			
				$tmpSums[$column-4]+=$actCov;	
				#$tmpCounts[$column-4]+=1;
			}
		
			if ($in_insertion){
				$lasts-=1;	
				if (!$lasts){$in_insertion=0}
			}	
			#print Dumper(\@tmpSums);
		}
		
		#print Dumper(\@tmpSums);
		

	}	
	 
} 
close fileHandle;

#print Dumper($totalLength);
#print Dumper(\@tmpSums);

my $cols = join(',', @POP_COLUMNS);
my $cols_print = join("\t", @POP_COLUMNS);

foreach my $column (@POP_COLUMNS){
	my $out= $tmpSums[$column-4]/$totalLength;
	print $out."\t";
}

print "\n";

=pod

=head1 NAME

coverage-for-sync-parts-prettybase.pl

=head1 SYNOPSIS

perl coverage-for-sync-parts-prettybase.pl --input syncPileupFileName --prettybase-summary fileName
--type 0/1/2 --populations columns 4,5,6

=head1 DESCRIPTION

The script calculates average coverage for regions defined in a prettybase-summary file. For the coverage calculation 
populational data from an input synchronized pileup file are used.

Output of the script is a single line that contains average coverage for each populational column specified by parameter --type. 
The output average coverages are in the same order as the columns specified by --type.  

=head1 OPTIONS

=over 4 

=item --input 

A synchronized pileup file. See script synchronize-pileup.pl in popoolation for details. The synchronized pileup 
file should contain only data for one chromosome. 

=item --prettybase-summary

A tab delimited file that contains following columns: position, 
type of feature ("INS" for insertion or "SNP"), significance ("SIG" or "NON"), length. You can create a prettybase-summary file from a prettybase file using script PopGenTools/misc/tina/summarize-prettybase.pl, see its help page for more details.

=item --type

A parameter that specifies whether to calculate average coverage for SNPs only (--type 0) or for insertions only (--type 1)
or for all together (--type 2 - default value of the parameter). 

=item --population

A comma separated list of columns thhat contain populations for which coverage is calculated. The columns are numbered from 1, 
first column containing population data is column 4 (see example)

=item --help, -h

Prints the help page.

=head1 EXAMPLE

=over 4

=item INPUT synchronized file example.sync 

	2L	1	G	0:1:15:90:0:0	0:0:25:100:0:0	0:0:80:20:0:0	0:0:95:15:0:0
	2L	6	C	0:1:15:90:0:0	0:0:25:100:0:0	0:0:80:20:0:0	0:0:95:15:0:0
	2L	10	C	0:0:80:1:0:0	0:0:115:0:0:0	0:1:80:3:0:0	0:1:95:0:0:0

=item INPUT prettybase-summary file example.summary

	1	SNP	SIG	1
	6	SNP	NON	1
	10	INS	NON	1

=item Used commands

perl coverage-for-sync-pileup-parts.pl --input example.sync --prettybase-summary example.summary --type 0 --population 4,5,6,7 

perl coverage-for-sync-pileup-parts.pl --input example.sync --prettybase-summary example.summary --type 1 --population 4,5,6,7 

perl coverage-for-sync-pileup-parts.pl --input example.sync --prettybase-summary example.summary --type 2 --population 4,5,6,7 

=item OUTPUTS

	106	125	100	110	
	81	115	84	96	
	97.6666666666667	121.666666666667	94.6666666666667	105.333333333333	

=back

=cut
