#!/usr/bin/perl-w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path;
#Date: 23-02-2010
# Author: Ram Vinay Pandey 

# Define the variables
my $inputfile = "";
my $outputfile = "";
my $number_of_species = 3;
my $alignment_dir = "";
my ($covered_length,$alignment_with_gap,$alignment_without_gap,$st,$nd) = (0,0,0,0,0);
my @startpos = ();
my @endpos = ();
my $help=0;
my $test=0;
my $verbose=1;



GetOptions(
    "alignment-dir=s"=>\$alignment_dir,
    "output=s"=>\$outputfile,
    "help" =>\$help,
    "test"=>\$test
) or pod2usage(-msg=>"Option(s) not valid $!",-verbose=>1);

pod2usage(-verbose=>2) if $help;
Test::runTests() if $test;

pod2usage(-msg=>"Provide Mauve alignment output directory path as input directory\n",-verbose=>1) unless $alignment_dir ;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $outputfile;
    


my ( $name, $path, $extension ) = File::Basename::fileparse ( $outputfile, '\..*' );
my $unaligned_region_output_file = $path."/".$name."-unaligned-regions.txt";

#rmtree($alignment_dir, 0, 0) if(-d $alignment_dir);

open my $ofh,">$outputfile" or die "Could not open file handle: $outputfile\n";
open my $uofh,">$unaligned_region_output_file" or die "Could not open file handle: $unaligned_region_output_file\n";


# Opening the mauve alignment files directory for reading.
opendir(DIR, $alignment_dir) || die ("Cannot open directory $alignment_dir in mauve parser $!\n");
my @files = readdir(DIR);
closedir(DIR);




print "\nMauve parsing started at\t". localtime() ." ....\n\n";

print "\nChromosome\tTotal aligned length\tAligned without gap (%age)\tAligned with gap (%age)\tUnaligned (%age)\n";
foreach my $file (@files) {
    
    unless ( ($file eq ".") || ($file eq "..") || ($file eq ".DS_Store")) {
	
	my $f = $file;
	my @f = split("-",$file);
	my $chromosome = $f[0];
	$inputfile = "$alignment_dir/$file";
	my ( $name, $path, $extension ) = File::Basename::fileparse ( $file, '\..*' );
	
	if ($extension eq ".alignment") {
	    
	    my $data = [];
	    ($alignment_with_gap,$alignment_without_gap,$st,$nd,$data) = MauveParser::Mauve_parser($inputfile,$outputfile,$chromosome,$ofh);
	    
	    my $chr_length = ($$nd[-1] - $$st[0])+1;
	    my $unaligned_region = MauveParser::calculate_unaligned_region($chromosome,$chr_length,$st,$nd,$uofh);
	    
	    #push (@$data,@$d);
	    my $stat = MauveParser::calculate_alignment_statistics($chr_length,$alignment_with_gap,$alignment_without_gap,$unaligned_region);
	    #my $stat = MauveParser::write_output($ofh,$data,$chr_length);
	    print "$chromosome\t$stat\n";
	    
	}
		
    }
    
}

print "\n\nMauve parsing completed at\t". localtime() ." ....\n";
exit;		


##########################################################################################
# MauveParser package contians many utility functions to parse the Mauve alignment output
# file and creates the tab delimeted output file.
##########################################################################################

{
    package MauveParser;
    
    sub write_output {
	
	my ($ofh,$data) = (shift,shift);
	my ($aligned_gap,$aligned) = (0,0);
	
	#[ sort {$a->{chr} cmp $b->{chr} or $a->{pos} <=> $b->{pos} } @$data];
	$data = [ sort { $a->{pos} <=> $b->{pos} } @$data];
	for my $d (@$data) {
	    
	    print $ofh "$d->{chr}\t$d->{pos}\t$d->{astate}\t$d->{dstate}\t$d->{sp1}\t$d->{sp2}\t$d->{sp3}\n";
	    
	    if (($d->{sp1} eq "-") or ($d->{sp2} eq "-") or ($d->{sp3} eq "-")) {
		$aligned_gap++;
	    }
	    else {
		$aligned++;
	    }
	    
	}

	return ($aligned,$aligned_gap);
    }
    

    sub calculate_alignment_statistics {
	
	my ($chr_length,$alignment_with_gap,$alignment_without_gap,$unaligned_region) = (shift,shift,shift,shift);
	
	my ($percent1,$percent2,$percent3) = (0,0,0);

        if ($alignment_without_gap>0) {
	    $percent1 = ($alignment_without_gap/$chr_length)*100;
	}
	
	if ($alignment_with_gap>0) {
	    $percent2 = ($alignment_with_gap/$chr_length)*100;
	}
	if ($unaligned_region>0) {
	    $percent3 = ($unaligned_region/$chr_length)*100;
	}
	
        $percent1 = sprintf "%.2f",$percent1;
        $percent2 = sprintf "%.2f",$percent2;
	$percent3 = sprintf "%.2f",$percent3;
	
	#print "$chr_length\t$aligned ($percent1)\t$aligned_gap ($percent2)\t$unaligned ($percent3)";
	
	my $stat = "";
	$stat = "$chr_length\t$alignment_without_gap ($percent1)\t$alignment_with_gap ($percent2)\t$unaligned_region ($percent3)";
	return $stat;
	
    }
    ##########################################################################################
    # Mauve_parser function parses the Mauve alignment output
    # file and creates the tab delimeted output file.
    ##########################################################################################
    sub Mauve_parser {
	
	my ($inputfile,$outputfile,$chr,$ofh) = (shift,shift,shift,shift);
	
	my $covered_length = 0;
	my $partner_count = 0;
	my @startpos = ();
	my @endpos = ();
	my $data = [];
	my ($alignment_with_gap,$alignment_without_gap) = (0,0);
	open (IN,"<$inputfile") || die ("could not open file $inputfile for reading $!\n");

	my $ct1 = 0;

	my @lines = ();
	while(my $line = <IN>) {
	    chomp $line;

	    # discard blank line
		if ($line =~ m/^\s*$/g) {
		    next;
		}
		# discard comment line
		elsif($line =~ m/^#/) {
		    next;
		}
		elsif ($line =~ m/=/) {
		    $ct1++;
		    foreach my $l(@lines) {
			if($l =~ m/\>/) {
			    $partner_count++;
			}
		    }

		    if($partner_count==$number_of_species) {
			
			my $aligned_block_sequence = MauveParser::remove_linebreak(\@lines,$outputfile);

			$covered_length += $$aligned_block_sequence[0]->{alignment_length};
			push (@startpos,$$aligned_block_sequence[0]->{startpos});
			push (@endpos,$$aligned_block_sequence[0]->{endpos});
			
			# Calculate the allelic states.
			my ($aligned,$aligned_gap) = MauveParser::calculate_allelic_state($aligned_block_sequence,$chr,$ofh);
			#my ($aligned,$aligned_gap) = MauveParser::write_output($ofh,$d);
			
			$alignment_with_gap += $aligned_gap;
			$alignment_without_gap += $aligned;
			
			#push (@$data,@$d);
		    }
		    
		    @lines = ();
		    $partner_count=0;
		}
		
		else {
		    push(@lines,$line);
		}
		
		#last if $ct1>1;
	}

	close IN;
	return ($alignment_with_gap,$alignment_without_gap,\@startpos,\@endpos);
    }
    
    
    ##########################################################################################
    # remove_linebreak function removed the newline and carriege return character from the
    # end of each line of Mauve alignment output file.
    ##########################################################################################
    sub remove_linebreak {
	
	my $alignment_block = shift;
	
	my @alignment_block = @{$alignment_block};
	my @alignment_block_sequences = ();
	push(@alignment_block,">");
	my @alignment_block_temp = ();


	foreach my $line(@alignment_block) {

	    if($line =~ m/\>/) {

		my $seq = "";
		if (scalar(@alignment_block_temp)>0) {
		    foreach my $l (@alignment_block_temp) {
			$seq .= $l;
		    }
		}

		push(@alignment_block_sequences,$seq);
		push(@alignment_block_sequences,$line);
		
		$seq = "";
		@alignment_block_temp = ();
		
	    }
	    else {
		push(@alignment_block_temp,"$line");
	    }
	}

	shift(@alignment_block_sequences);
	pop(@alignment_block_sequences);
	
	my $ct=0;

	my $header = "";

	my $aligned_block_sequence = [];
	
	foreach my $line (@alignment_block_sequences) {
	    chomp $line;
	    # discard blank line
	    if ($line =~ m/^\s*$/g) {
		next;
	    }
	    # discard comment line
	    elsif($line =~ m/^#/) {
		next;
	    }
	    elsif($line =~ m/\>/) {
		$header = $line;	    
	    }
	    else {
		#> 1:26790-236196 + /Volumes/Main/fst-project-05-01-2010/fisher-exact-test/Mauve-analysis/chr-input-flybase/4-dmel.fasta
		my @fasta_header = split(" ",$header);
		
		my $strand = $fasta_header[2];
		
		$fasta_header[1] =~ m/(.*)\:(.*)\-(.*)/;
		my ($start,$end,$length) = ($2,$3,(abs($3-$2))+1);
		
		my ( $filename, $path, $extension ) = File::Basename::fileparse ( $fasta_header[3], '\..*' );
		push @$aligned_block_sequence,
		{
		seqheader=>$header,
		sequence=>$line,
		strand=>$strand,
		file=>$filename,
		startpos=>$start,
		endpos=>$end,
		alignment_length=>$length
		};
	    }
	}

	return $aligned_block_sequence;
    }
    
    
    ##########################################################################################
    # calculate_allelic_state function calculates the allelic state in species 1, 2 and 3
    # at each locus within alignment block.
    ##########################################################################################
    
    sub calculate_allelic_state {
	my ($aligned_block_sequence,$chr) = (shift,shift,$ofh);
	my $data=[];
	my $alleles = [];
	
	my @seq1 = split("",$$aligned_block_sequence[0]->{sequence});
	my @seq2 = split("",$$aligned_block_sequence[1]->{sequence});
	my @seq3 = split("",$$aligned_block_sequence[2]->{sequence});
	my ($aligned_gap,$aligned) = (0,0);
	
	my $startpos = $$aligned_block_sequence[0]->{startpos};

	for(my $i=0; $i<scalar(@seq1);$i++) {
	    
	    if (($seq1[$i] eq "-") or ($seq2[$i] eq "-") or ($seq3[$i] eq "-")) {
		$aligned_gap++;
	    }
	    else {
		$aligned++;
	    }
		
	    if (($seq1[$i] ne "-") and ($seq2[$i] ne "-") and ($seq3[$i] ne "-")) {
		push @$alleles,
		my $alleles = {
		    sp1=>$seq1[$i],
		    sp2=>$seq2[$i],
		    sp3=>$seq3[$i]
		};
		
		my ($ancestral_state,$derived_state) = MauveParser::check_allelic_state($alleles);
		print $ofh "$chr\t$startpos\t$ancestral_state\t$derived_state\t$seq1[$i]\t$seq2[$i]\t$seq3[$i]\n";
		
		
		$alleles = {};
		
		#print $ofh "$chr\t$startpos\t$ancestral_state\t$derived_state\t$seq1[$i]\t$seq2[$i]\t$seq3[$i]\n";
		$startpos++;
	    }
	}
	return ($aligned,$aligned_gap);
    } # Closed function calculate_allelic_state
    
    
    ##########################################################################################
    # check_allelic_state function calculates the parental and derived allelic state
    # based upon 3 species allelic state at each locus within alignment block.
    ##########################################################################################
    
    sub check_allelic_state {
	
	my $alleles = shift;
	# $alleles->{sp1}    = Drosophila melanogaster
	# $alleles->{sp2}   =  Drosophila simulans
	# $alleles->{sp3}   =  Drosophila yakuba
	
	my ($ancestral_state,$derived_state) = ("","");
	
	if(($alleles->{sp1} =~ m/N/gi) or ($alleles->{sp2} =~ m/N/gi) or ($alleles->{sp3} =~ m/N/gi)) {
	    ($ancestral_state,$derived_state) = ("na","na");
	}
	elsif(($alleles->{sp1} eq "-") or ($alleles->{sp2} eq "-") or ($alleles->{sp3} eq "-")) {
	    ($ancestral_state,$derived_state) = ("na","na");
	}
	else {
	    if (("$alleles->{sp2}" eq "$alleles->{sp3}") and ("$alleles->{sp1}" ne "$alleles->{sp2}")) {
		($ancestral_state,$derived_state) = ("$alleles->{sp3}","$alleles->{sp1}");
	    }
	    elsif (("$alleles->{sp1}" eq "$alleles->{sp3}") and ("$alleles->{sp1}" ne "$alleles->{sp2}")) {
		($ancestral_state,$derived_state) = ("$alleles->{sp1}","0");
	    }
	    elsif (("$alleles->{sp2}" eq "$alleles->{sp3}") and ("$alleles->{sp1}" eq "$alleles->{sp2}")) {
		($ancestral_state,$derived_state) = ("$alleles->{sp1}","0");
	    }
	    else {
		($ancestral_state,$derived_state) = ("na","na");
	    }
	}
	return ($ancestral_state,$derived_state);
    }
    
    
    
    ##########################################################################################
    # calculate_unaligned_region function calculates the unaligned region with respect to
    # reference chromosome.
    ##########################################################################################
    
    sub calculate_unaligned_region {
	my ($chromosome,$chr_length,$st,$nd,$ofh) = (shift,shift,shift,shift,shift);
	my @startpos = @{$st};
	my @endpos = @{$nd};
	my $data = [];
	# Calculating the unaligned region of given chromosome

	shift(@startpos);
	my $last_endpos = pop(@endpos);
	
	my $unaligned_region = 0;
	for (my $i=0;$i<scalar(@startpos);$i++) {
	    my $st1 = $endpos[$i]+1;
	    my $nd1 = $startpos[$i]-1;
	    
	    my $length = ($nd1-$st1)+1;
	    $unaligned_region += $length;
	    
	    my $lengthKB = $length/1000;
	    $lengthKB = sprintf "%.2f",$lengthKB;
	    print $ofh "$chromosome\t$st1\t$nd1\t$length\t$lengthKB\n";
=cut
	    for my $i ($st1 .. $nd1) {
	    push @$data,
		{
		    chr=>$chromosome,
		    pos=>$i,
		    astate=>0,
		    dstate=>0,
		    sp1=>0,
		    sp2=>0,
		    sp3=>0
		};
	    }
=cut
	    
	}
	return $unaligned_region;
	
    }
    

} # MauveParser package completed


#
#    "chralignment-dir=s"=>\$alignment_dir,
#    "output=s"=>\$outputfile,
#    "help" =>\$help,
#    "test"=>\$test
    
=head1 NAME

mauve-parser.pl - Parses the Mauve chromosome alignement output file for given number of chromosomes.

=head1 SYNOPSIS

 perl mauve-parser.pl --alignment-dir mauve-alignment-files-dir --output mauve-parser-output.txt

=head1 OPTIONS

=over 4

=item B<--alignment-dir>

the input directory/folder full path where Mavue alignment output for one or multiple chromosome are stored: "file name should be Chromosome-alignment.alignment", example: 2L-alignment.alignment.


=item B<--output>

The output file. Mandatory parameter


=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

Input is a folder/directoy which contains mauve alignment output for all chromosomes. Within input directory each mauve multiple alignment output file looks like following:

 #FormatVersion Mauve1
 #Sequence1File	2L-dmel.fasta
 #Sequence1Format	FastA
 #Sequence2File	Kib32_sorted_merged_2L.fa
 #Sequence2Format	FastA
 #Sequence3File	2L-dyak.fasta
 =
 > 1:1781127-1781186 + 2L-dmel.fasta
 GCTGTCAGCAATTTGATTTGATTTGCCAACTCTGCAAAAATATTTCATAAATATTTGGCC
 > 2:1762817-1762876 - Kib32_sorted_merged_2L.fa
 GCTGTCAGCAATTTGATTTGATTTGCCAACTCTGCAAAAATATTTCATAAATATTTGCCC
 > 3:1756489-1756548 + 2L-dyak.fasta
 GCTGTCAGCAATTTGATTTGATTTGCCAACTCTGCGAAAATATTTCATAAATATTTGCCC
 =
 First sequence should be of reference chromosome (Example D. melanogaster chromosome sequence)
 Second sequence should be of intermidate species chromosome (Example D. simulans chromosome sequence)
 Third sequence should be of outgroup species chromosome (Example D. yakuba chromosome sequence)
 
=head2 Output

An output of this program looks like in the given example:

 2L	5783	na	na	G	G	G
 2L	5784	na	na	C	C	C
 2L	5785	na	na	C	T	G

 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: ancestral allelic state
 col 4: derived allelic state
 col 5: allelic state in species 1 (reference species)
 col 6: allelic state in species 2
 col 7: allelic state in species 3 (out group species)

=head1 AUTHORS

Ram vinay pandey

Robert Kofler

Pablo Orozco terWengel

Christian Schloetterer

=cut