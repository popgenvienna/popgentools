#!/usr/bin/perl

package HkaTestModule;

use warnings;
use strict;

use FindBin qw($Bin);
#use lib "$Bin/External/Statistics/"; 
use lib "$Bin/External/";

#PopGenTools
use Pileup;
use GffGtfParser;
use VarianceExactCorrection;
use VarMath qw(get_thetadiv_buffer get_thetadiv_buffer_sqr);

use Data::Dumper;
use Statistics::Distributions qw(chisqrprob);

#export
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
);
#upgrade to more than one param set?

# not tested/used functionality:
#	
#
#
#
our @EXPORT_OK = qw(
	array_to_tab_separ_string
	calculate_and_store_chi_square_statistic_for_each_gene
	calculate_and_store_chi_square_statistic_for_one_gene
	calculate_and_store_D_expectations_and_variance
	calculate_and_store_measures_for_genes
	calculate_and_store_T
	contains
	create_actual_list_of_genes_on_position
	create_list_of_genes_with_non_zero_S_variance
	create_list_of_long_enough_genes
	create_list_of_nonoverlapping_genes
	find_overlaps_for_genes
	handling_with_overlapping_and_short_and_zero_variance_genes
	load_all_data_GTF_PILEUP_MAUVEparsed
	load_genes_measures
	load_GPM_collect_genes_calculate_measures_and_overlaps
	load_GPM_collect_genes_calculate_measures_and_overlaps_old_theta
	load_MAUVEparsed_2_reference_seqs_well_defined_annotated
	load_overlapping_genes
	load_PILEUP_coveredPositions_annotated_mauved
	load_previously_calculated_measures_and_overlaps
	load_reference_gene_list
	load_simulated_GP_collect_genes_calculate_measures_and_overlaps
	load_simulated_GP_collect_genes_calculate_measures_and_overlaps_old_theta
	max
	print_first_round_of_data
	print_genes_measures
	print_log_multiple_geneIDs_positions
	print_log_removed_genes
	print_log_splitted_genes
	print_log_withdrawn
	print_overlapping_genes
	print_program_params_info
	print_reference_set_of_genes_with_T_for_the_set
	print_selected_genes_with_all_statistics
	print_selected_genes_without_T
	randomly_choose_genes
	select_CDS_gene_well_defined_parts
	select_the_effectively_longest_transcript_for_each_gene
	set_strand_for_gene
	store_genes_info
	update_length_of_gene
	update_list_of_chromosomes
	update_measure_sums
	update_min_max_bounds_of_gene
	update_strand_values
);

#########################
##					   ##	
## some misc functions ##
##                     ##
#########################

# subroutine max returns the max value from the input list of numbers
sub max{
	my $m=$_[0];
	for (my $i=0; $i< scalar(@_); $i++){
		next unless ($_[$i]>$m);
		$m = $_[$i];		
	} 
	return $m;
}

# subroutine contains returns bool:
# true if $ptrArray contains an element $geneID, false otherwise
sub contains{
	my ($geneID, $ptrArray) = @_;
	my $is_in_array=0;
		
	foreach my $element (@$ptrArray){
		next unless ($element eq $geneID);
		$is_in_array = 1;		
	}
	return $is_in_array;
}


#######################################################
####    										   ####	
#### loading data from GTF, PILEUP and MAUVE files ####
####											   ####				
#######################################################

sub load_all_data_GTF_PILEUP_MAUVEparsed{
	my ($GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
	$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL
	)=@_;
	
	my $ptrGeneIDsAllData = {};
	my $ptrGeneParts = {};
	
	my $ptrAnnotation = {};
	my $ptrReferences = {};
	my $ptrCoveredPositions = {};
	my $ptrWithdrawnInGtf = {};
	
	
	($ptrAnnotation, $ptrWithdrawnInGtf) = load_GTF_file_CDS_with_gene_IDs_skip_withdrawn_genes($GTF_FILE, $WITHDRAWN_FILE, $ptrGeneIDsAllData);
	($ptrAnnotation, $ptrReferences) = load_MAUVEparsed_2_reference_seqs_well_defined_annotated($MAUVE_PARSED_FILE, $ptrAnnotation);
	($ptrAnnotation, $ptrReferences, $ptrCoveredPositions) = load_PILEUP_coveredPositions_annotated_mauved($PILEUP_FILE, $QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL, $ptrAnnotation, $ptrReferences);
	
	my $ptrAllData={};
	
	foreach my $chromosome (keys %{$ptrAnnotation}){
		my $end=scalar(@{$ptrAnnotation->{$chromosome}})-1;

		for (my $i=0; $i<$end; $i++){
			
			#print Dumper($ptrAnnotation, $ptrReferences, $ptrCoveredPositions);
			
			next unless ( defined($ptrAnnotation->{$chromosome}[$i]) );
			#next unless ( (defined($ptrAnnotation->{$chromosome}[$i])) and (exists($ptrCoveredPositions->{$chromosome}[$i])));
			#next unless (exists($ptrCoveredPositions->{$chromosome}[$i]));
			
			$ptrAllData->{$chromosome}[$i]{feature}=$ptrAnnotation->{$chromosome}[$i]{feat};
			$ptrAllData->{$chromosome}[$i]{strand}= $ptrAnnotation->{$chromosome}[$i]{strand};
			$ptrAllData->{$chromosome}[$i]{geneID}=$ptrAnnotation->{$chromosome}[$i]{geneID};
			#store codon position;
			$ptrAllData->{$chromosome}[$i]{codonPosition}=$ptrAnnotation->{$chromosome}[$i]{codonPosition};
			
			$ptrAllData->{$chromosome}[$i]{secondRef}= $ptrReferences->{$chromosome}[$i]{secondRef};
			
			$ptrAllData->{$chromosome}[$i]{firstRef} = $ptrCoveredPositions->{$chromosome}[$i]{firstRef};
			$ptrAllData->{$chromosome}[$i]{allele1} = $ptrCoveredPositions->{$chromosome}[$i]{allele1};
			#$ptrAllData->{$chromosome}[$position]{allele1count} = $ptrCoveredPositions->{$chromosome}[$position]{allele1count};
			$ptrAllData->{$chromosome}[$i]{allele2} = $ptrCoveredPositions->{$chromosome}[$i]{allele2};
			#$ptrAllData->{$chromosome}[$position]{allele2count} = $ptrCoveredPositions->{$chromosome}[$position]{allele2count}; 		
			
			$ptrAllData->{$chromosome}[$i]{eucov}=$ptrCoveredPositions->{$chromosome}[$i]{eucov};
			$ptrAllData->{$chromosome}[$i]{ispuresnp}=$ptrCoveredPositions->{$chromosome}[$i]{ispuresnp};
			$ptrAllData->{$chromosome}[$i]{A}=$ptrCoveredPositions->{$chromosome}[$i]{A};
			$ptrAllData->{$chromosome}[$i]{T}=$ptrCoveredPositions->{$chromosome}[$i]{T};
			$ptrAllData->{$chromosome}[$i]{C}=$ptrCoveredPositions->{$chromosome}[$i]{C};
			$ptrAllData->{$chromosome}[$i]{G}=$ptrCoveredPositions->{$chromosome}[$i]{G};
			$ptrAllData->{$chromosome}[$i]{chr}=$ptrCoveredPositions->{$chromosome}[$i]{chr};
			$ptrAllData->{$chromosome}[$i]{pos}=$ptrCoveredPositions->{$chromosome}[$i]{pos};
			
						
		}
	}
	$ptrGeneParts = select_CDS_gene_well_defined_parts($ptrAllData, $ptrGeneIDsAllData);	
	return ($ptrAllData, $ptrGeneParts, $ptrWithdrawnInGtf);
}


#sub load_all_data_GTF_PILEUP_MAUVEparsed_old_theta{
#	my ($GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
#	$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL
#	)=@_;
#	
#	my $ptrGeneIDsAllData = {};
#	my $ptrGeneParts = {};
#	
#	my $ptrAnnotation = {};
#	my $ptrReferences = {};
#	my $ptrCoveredPositions = {};
#	my $ptrWithdrawnInGtf = {};
#	
#	($ptrAnnotation, $ptrWithdrawnInGtf) = load_GTF_file_CDS_with_gene_IDs_skip_withdrawn_genes($GTF_FILE, $WITHDRAWN_FILE, $ptrGeneIDsAllData);
#	($ptrAnnotation, $ptrReferences) = load_MAUVEparsed_2_reference_seqs_well_defined_annotated($MAUVE_PARSED_FILE, $ptrAnnotation);
#	($ptrAnnotation, $ptrReferences, $ptrCoveredPositions) = load_PILEUP_coveredPositions_annotated_mauved_old_theta($PILEUP_FILE, $QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL, $ptrAnnotation, $ptrReferences);
#	
#	my $ptrAllData={};
#	
#	foreach my $chromosome (keys %{$ptrAnnotation}){
#		my $end=scalar(@{$ptrAnnotation->{$chromosome}})-1;
#
#		for (my $i=0; $i<$end; $i++){
#			
#			#print Dumper($ptrAnnotation, $ptrReferences, $ptrCoveredPositions);
#			
#			next unless ( defined($ptrAnnotation->{$chromosome}[$i]) );
#			#next unless ( (defined($ptrAnnotation->{$chromosome}[$i])) and (exists($ptrCoveredPositions->{$chromosome}[$i])));
#			#next unless (exists($ptrCoveredPositions->{$chromosome}[$i]));
#			
#			$ptrAllData->{$chromosome}[$i]{feature}=$ptrAnnotation->{$chromosome}[$i]{feat};
#			$ptrAllData->{$chromosome}[$i]{strand}=$ptrAnnotation->{$chromosome}[$i]{strand};
#			$ptrAllData->{$chromosome}[$i]{geneID}=$ptrAnnotation->{$chromosome}[$i]{geneID};
#			$ptrAllData->{$chromosome}[$i]{codonPosition}=""
#			
#			$ptrAllData->{$chromosome}[$i]{offset}=$ptrAnnotation->{$chromosome}[$i]{offset};
#			#$ptrAllData->{$chromosome}[$i]{firstRef}= $ptrReferences->{$chromosome}[$i]{firstRef};
#			$ptrAllData->{$chromosome}[$i]{secondRef}= $ptrReferences->{$chromosome}[$i]{secondRef};
#			$ptrAllData->{$chromosome}[$i]{ispuresnp}=$ptrCoveredPositions->{$chromosome}[$i]{ispuresnp};
#			$ptrAllData->{$chromosome}[$i]{eucov}=$ptrCoveredPositions->{$chromosome}[$i]{eucov};
#			$ptrAllData->{$chromosome}[$i]{A}=$ptrCoveredPositions->{$chromosome}[$i]{A};
#			$ptrAllData->{$chromosome}[$i]{T}=$ptrCoveredPositions->{$chromosome}[$i]{T};
#			$ptrAllData->{$chromosome}[$i]{C}=$ptrCoveredPositions->{$chromosome}[$i]{C};
#			$ptrAllData->{$chromosome}[$i]{G}=$ptrCoveredPositions->{$chromosome}[$i]{G};
#			$ptrAllData->{$chromosome}[$i]{chr}=$ptrCoveredPositions->{$chromosome}[$i]{chr};
#			$ptrAllData->{$chromosome}[$i]{pos}=$ptrCoveredPositions->{$chromosome}[$i]{pos};
#			
#		}
#	}
#	$ptrGeneParts = select_CDS_gene_well_defined_parts($ptrAllData, $ptrGeneIDsAllData);	
#	return ($ptrAllData, $ptrGeneParts, $ptrWithdrawnInGtf);
#}



sub load_MAUVEparsed_2_reference_seqs_well_defined_annotated{
	my ($IN_FILE, $ptrAnnotation)=@_;
	
	my $ptrReferences={};
	
	open MAUVE, "<", $IN_FILE or die "Could not open mauve-parsed file $IN_FILE";
	while (my $line = <MAUVE>){
		chomp($line);
		my ($chromosome, $position, $first, $second) = split "\t", $line;
 		
 		next unless defined($ptrAnnotation->{$chromosome}[$position]);
 		
 		if ( $second =~ m/[ATCG]/ ){
 			#$ptrReferences->{$chromosome}[$position]{firstRef}=$first;
			$ptrReferences->{$chromosome}[$position]{secondRef}=$second;
 		}else{
 			$ptrAnnotation->{$chromosome}[$position]=undef;
 			$ptrReferences->{$chromosome}[$position]=undef;
 		}
		
	}
	close MAUVE;
	return ($ptrAnnotation, $ptrReferences);
}


sub load_PILEUP_coveredPositions_annotated{
	my ($IN_FILE,
		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
		$ptrAnnotation
	) = @_;
	
	my $pileupParser=get_pileup_parser($QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL);
	#my $pileupParserOld=get_pileup_parser($QUAL_ENCODING, 1, $MIN_COV, $MAX_COV, $MIN_QUAL);
	my $ptrCovered={};
	
	open PILEUP, "<", $IN_FILE or die "Could not open pileup file $IN_FILE.";
	while (my $line = <PILEUP>){
		chomp($line);
		
		my $parsedLine = $pileupParser->($line);
		#my $parsedLineOld = $pileupParserOld->($line);
		my $chromosome=$parsedLine->{chr};
		my $position=$parsedLine->{pos};
		
		next unless (defined($ptrAnnotation->{$chromosome}[$position]));
		next unless ($parsedLine->{iscov}); 
		$ptrCovered->{$chromosome}[$position]=$parsedLine;

		
	}	
	close PILEUP;	
	
	foreach my $chr (keys %$ptrAnnotation){ 
		for (my $i=0; $i< scalar(@{$ptrAnnotation->{$chr}}); $i++){
			next if (exists($ptrCovered->{$chr}[$i]));
			$ptrCovered->{$chr}[$i] = undef;
			$ptrAnnotation->{$chr}[$i]=undef;
		}
	}
	return $ptrCovered;
}

sub load_PILEUP_coveredPositions_annotated_mauved{
	my ($IN_FILE,
		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
		$ptrAnnotation, $ptrReferences
	) = @_;
	
	my $pileupParser=get_pileup_parser($QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL);
	#my $pileupParserOld=get_pileup_parser($QUAL_ENCODING, 1, $MIN_COV, $MAX_COV, $MIN_QUAL);
	my $ptrCovered={};
	
	open PILEUP, "<", $IN_FILE or die "Could not open pileup file $IN_FILE.";
	while (my $line = <PILEUP>){
		chomp($line);
		
		my $parsedLine = $pileupParser->($line);
		#my $parsedLineOld = $pileupParserOld->($line);
		my $chromosome=$parsedLine->{chr};
		my $position=$parsedLine->{pos};
		
		next unless ( (defined($ptrAnnotation->{$chromosome}[$position])) and (defined($ptrReferences->{$chromosome}[$position])) );
		next unless ($parsedLine->{iscov}); 
		$ptrCovered->{$chromosome}[$position]=$parsedLine;
#choose major allele
		my %tmpCounts = ();
		
		$tmpCounts{"A"} = $parsedLine->{"A"};
		$tmpCounts{"T"} = $parsedLine->{"T"};
		$tmpCounts{"C"} = $parsedLine->{"C"};
		$tmpCounts{"G"} = $parsedLine->{"G"};
		
		my $ptrSortedByValue=[];
		
		foreach my $nucleotide (sort {$tmpCounts{$b}<=>$tmpCounts{$a}} (keys(%tmpCounts))){
			my %hash = (nucleotide=>$nucleotide, count=>$tmpCounts{$nucleotide});
			push @{$ptrSortedByValue}, \%hash;
		}
		
		$ptrCovered->{$chromosome}[$position]{firstRef} = $ptrSortedByValue->[0]{nucleotide};
		
		if ($parsedLine->{ispuresnp}){
			$ptrCovered->{$chromosome}[$position]{allele1} = $ptrSortedByValue->[0]{nucleotide};
			$ptrCovered->{$chromosome}[$position]{allele2} = $ptrSortedByValue->[1]{nucleotide};
		} 
	}	
	close PILEUP;	
	
	foreach my $chr (keys %$ptrAnnotation){ 
		for (my $i=0; $i< scalar(@{$ptrAnnotation->{$chr}}); $i++){
			next if (exists($ptrCovered->{$chr}[$i]));
			$ptrCovered->{$chr}[$i] = undef;
			$ptrAnnotation->{$chr}[$i]=undef;
			$ptrReferences->{$chr}[$i]=undef;
		}
	}
	return ($ptrAnnotation, $ptrReferences, $ptrCovered);
}

#sub load_PILEUP_coveredPositions_annotated_mauved_old_theta{
#	my ($IN_FILE,
#		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
#		$ptrAnnotation, $ptrReferences
#	) = @_;
#	
#	my $pileupParser=get_pileup_parser($QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL);
#	my $ptrCovered={};
#	
#	open PILEUP, "<", $IN_FILE or die "Could not open pileup file $IN_FILE.";
#	while (my $line = <PILEUP>){
#		chomp($line);
#		
#		my $parsedLine = $pileupParser->($line);
#		my $chromosome=$parsedLine->{chr};
#		my $position=$parsedLine->{pos};
#		
#		next unless ( (defined($ptrAnnotation->{$chromosome}[$position])) and (defined($ptrReferences->{$chromosome}[$position])) );
#		next unless ($parsedLine->{iscov}); 
#		$ptrCovered->{$chromosome}[$position]=$parsedLine;
#		
#	}	
#	close PILEUP;	
#	
#	foreach my $chr (keys %$ptrAnnotation){ 
#		for (my $i=0; $i< scalar(@{$ptrAnnotation->{$chr}}); $i++){
#			next if (exists($ptrCovered->{$chr}[$i]));
#			$ptrCovered->{$chr}[$i] = undef;
#			$ptrAnnotation->{$chr}[$i]=undef;
#			$ptrReferences->{$chr}[$i]=undef;
#		}
#	}
#	return ($ptrAnnotation, $ptrReferences, $ptrCovered);
#}



sub load_simulated_all_data_GTF_PILEUP{
		my ($GTF_FILE, $PILEUP_FILE, 
	$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL
	)=@_;
	
	
	my $ptrGeneIDsAllData = {};
	my $ptrGeneParts = {};
	
	my ($ptrAnnotation) = load_GTF_file_CDS_with_gene_IDs($GTF_FILE, $ptrGeneIDsAllData);
	my $ptrCoveredPositions = load_PILEUP_coveredPositions_annotated($PILEUP_FILE, $QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL, $ptrAnnotation);
	
	my $ptrAllData={};
	
	foreach my $chromosome (keys %{$ptrAnnotation}){
		my $end=scalar(@{$ptrAnnotation->{$chromosome}})-1;

		for (my $i=0; $i<$end; $i++){
			
			#print Dumper($ptrAnnotation, $ptrReferences, $ptrCoveredPositions);
			
			next unless ( defined($ptrAnnotation->{$chromosome}[$i]) );
			#next unless ( (defined($ptrAnnotation->{$chromosome}[$i])) and (exists($ptrCoveredPositions->{$chromosome}[$i])));
			#next unless (exists($ptrCoveredPositions->{$chromosome}[$i]));
			
			$ptrAllData->{$chromosome}[$i]{feature}=$ptrAnnotation->{$chromosome}[$i]{feat};
			$ptrAllData->{$chromosome}[$i]{strand}=$ptrAnnotation->{$chromosome}[$i]{strand};
			$ptrAllData->{$chromosome}[$i]{geneID}=$ptrAnnotation->{$chromosome}[$i]{geneID};
			$ptrAllData->{$chromosome}[$i]{eucov}=$ptrCoveredPositions->{$chromosome}[$i]{eucov};
			$ptrAllData->{$chromosome}[$i]{ispuresnp}=$ptrCoveredPositions->{$chromosome}[$i]{ispuresnp};
			$ptrAllData->{$chromosome}[$i]{firstRef}= "A";
			$ptrAllData->{$chromosome}[$i]{secondRef}= "A";
			$ptrAllData->{$chromosome}[$i]{A}=$ptrCoveredPositions->{$chromosome}[$i]{A};
			$ptrAllData->{$chromosome}[$i]{T}=$ptrCoveredPositions->{$chromosome}[$i]{T};
			$ptrAllData->{$chromosome}[$i]{C}=$ptrCoveredPositions->{$chromosome}[$i]{C};
			$ptrAllData->{$chromosome}[$i]{G}=$ptrCoveredPositions->{$chromosome}[$i]{G};
			$ptrAllData->{$chromosome}[$i]{chr}=$ptrCoveredPositions->{$chromosome}[$i]{chr};
			$ptrAllData->{$chromosome}[$i]{pos}=$ptrCoveredPositions->{$chromosome}[$i]{pos};
			
		}
	}
	$ptrGeneParts = select_CDS_gene_well_defined_parts($ptrAllData, $ptrGeneIDsAllData);	
	return ($ptrAllData, $ptrGeneParts);
}


sub load_simulated_GP_collect_genes_calculate_measures_and_overlaps{
	my (
		$GTF_FILE, $PILEUP_FILE, 
		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
		$POOL_SIZE, 
		$SEGREGATING_SITES
	)=@_;
		
	my ($ptrAllData, $ptrGeneParts) = load_simulated_all_data_GTF_PILEUP(
											$GTF_FILE, $PILEUP_FILE, 
											$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL
										);
	
	my ($ptrGeneData, $ptrSplittedGenes) = store_genes_info($ptrGeneParts, $ptrAllData);		
	my $ptrGeneIDs = select_the_effectively_longest_transcript_for_each_gene($ptrGeneData);
	
	my $ptrOverlaps={}; my $ptrOverlappingGeneIDs={};
	($ptrOverlaps, $ptrOverlappingGeneIDs, $ptrGeneIDs) = find_overlaps_for_genes($ptrAllData, $ptrGeneIDs);

	calculate_and_store_measures_for_genes($ptrGeneData, $ptrAllData, $ptrGeneIDs, $ptrGeneParts, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV, $SEGREGATING_SITES);
	calculate_and_store_S_exp_and_var_for_genes($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData, $POOL_SIZE, $MIN_COUNT, $SEGREGATING_SITES);
	
	#($ptrGeneData, $ptrGeneIDs, $POOL_SIZE, $MIN_COUNT)
	return ($ptrGeneIDs, $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs);		
}


#sub load_simulated_GP_collect_genes_calculate_measures_and_overlaps_old_theta{
#	my (
#		$GTF_FILE, $PILEUP_FILE, 
#		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
#		$POOL_SIZE, 
#		$SEGREGATING_SITES
#	)=@_;
#		
#	my ($ptrAllData, $ptrGeneParts) = load_simulated_all_data_GTF_PILEUP(
#											$GTF_FILE, $PILEUP_FILE, 
#											$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL
#										);
#	
#	my ($ptrGeneData, $ptrSplittedGenes) = store_genes_info($ptrGeneParts, $ptrAllData);		
#	my $ptrGeneIDs = select_the_effectively_longest_transcript_for_each_gene($ptrGeneData);
#	
#	my $ptrOverlaps={}; my $ptrOverlappingGeneIDs={};
#	($ptrOverlaps, $ptrOverlappingGeneIDs, $ptrGeneIDs) = find_overlaps_for_genes($ptrAllData, $ptrGeneIDs);
#
#	calculate_and_store_measures_for_genes_old_theta($ptrGeneData, $ptrAllData, $ptrGeneIDs, $ptrGeneParts, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV);
#	calculate_and_store_S_exp_and_var_for_genes($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData, $POOL_SIZE, $MIN_COUNT, $SEGREGATING_SITES);
#	return ($ptrGeneIDs, $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs);		
#}



sub set_syn_nonsyn_nothing_for_position{
	my ($ptrAllData, $ptrPositionData, $chromosome, $position, $start, $end, $strand, $ptrCodons, $alphabet)=@_;
	
	my $codonPosition = $ptrPositionData->{codonPosition};
	
	my $is_valid_seg_position = undef;
	#check position validity
	#all triplet positions well defined in the same part -- compare position to start, end
	if ($strand eq "+"){
		#good situations
		if (($codonPosition == 0) && ($position+2<=$end)){$is_valid_seg_position = 1}
		elsif (($codonPosition == 1) && ($position+1<=$end) && ($position-1>=$start)){ $is_valid_seg_position = 1}
		elsif (($codonPosition == 2) && ($position-2>=$start)){$is_valid_seg_position = 1}
		else{$is_valid_seg_position = 0}
	}elsif($strand eq "-"){
		if (($codonPosition == 0) && ($position-2>=$start)){$is_valid_seg_position = 1}
		elsif (($codonPosition == 1) && ($position-1>=$start) && ($position+1<=$end)){$is_valid_seg_position = 1}
		elsif (($codonPosition == 2) && ($position+2<=$end)){$is_valid_seg_position = 1}
		else{$is_valid_seg_position = 0} 	
	}

	my $allele1 = $ptrPositionData->{allele1};
	my $allele2 = $ptrPositionData->{allele2};
	
	if (!(defined($allele1) and defined($allele2))){
		$is_valid_seg_position = 0;
	}

	if ($is_valid_seg_position){
		my $triplet1;
		my $triplet2;
		my $aa1;
		my $aa2;

		my $ref_min_2 = $ptrAllData->{$chromosome}[$position-2]{firstRef};
		my $ref_min_1 = $ptrAllData->{$chromosome}[$position-1]{firstRef};
		my $ref_plus_1 = $ptrAllData->{$chromosome}[$position+1]{firstRef};
		my $ref_plus_2 = $ptrAllData->{$chromosome}[$position+2]{firstRef};

		if ($strand eq "+"){
			#translate triplets with 2 majors alleles at the position to amino acids
			#if the amino acids are the same, set ptrAllData->{$chromosome}[$position]{syn} = 1	
			if ($codonPosition == 0){
				$triplet1 = $allele1.$ref_plus_1.$ref_plus_2;
				$triplet2 = $allele2.$ref_plus_1.$ref_plus_2;
			}elsif($codonPosition == 1){
				$triplet1 = $ref_min_1.$allele1.$ref_plus_1;
				$triplet2 = $ref_min_1.$allele2.$ref_plus_1;
			}elsif($codonPosition == 2){
				$triplet1 = $ref_min_2.$ref_min_1.$allele1;
				$triplet2 = $ref_min_2.$ref_min_1.$allele2;				
			}			

		}elsif($strand eq "-"){
			#reverse-complement if strand eq "-"  
			#translate triplets with 2 majors alleles at the position to amino acids
			#if the amino acids are the same, set ptrAllData->{$chromosome}[$position]{syn} = 1
			if ($codonPosition == 0){
				$triplet1 = translate($allele1).translate($ref_min_1).translate($ref_min_2);
				$triplet2 = translate($allele2).translate($ref_min_1).translate($ref_min_2);
			}elsif($codonPosition == 1){
				$triplet1 = translate($ref_plus_1).translate($allele1).translate($ref_min_1);
				$triplet2 = translate($ref_plus_1).translate($allele2).translate($ref_min_1);				
			}elsif($codonPosition == 2){
				$triplet1 = translate($ref_plus_2).translate($ref_plus_1).translate($allele1);
				$triplet2 = translate($ref_plus_2).translate($ref_plus_1).translate($allele2);
			}
		}else{
			die "problem with strands\n";
		}

		#print Dumper($triplet1, $triplet2, $aa1, $aa2, $ptrCodons, join("", values %$ptrCodons));
		
		
		if (($triplet1 =~ m/[ATCG]{3}/) and ($triplet2 =~ m/[ATCG]{3}/)){
			$aa1 = $ptrCodons->{$triplet1};
			$aa2 = $ptrCodons->{$triplet2};
			
			#my $AAalphabet = "ILVFMCAGPTSYWQNHEDKR-";
			#if  ( ($aa1 =~ m/[\Qjoin("", values %$ptrCodons)\E]/) && ($aa2 =~ m/[\Qjoin("", values %$ptrCodons)\E]/)) {
			if  ( ($aa1 =~ m/[\Q$alphabet\E]/) && ($aa2 =~ m/[\Q$alphabet\E]/)) {
				if ($aa1 eq $aa2){	$ptrAllData->{$chromosome}[$position]{syn}=1	}
				else{	 $ptrAllData->{$chromosome}[$position]{syn}=0	}
			}else{
				die "problem with coding of amino acids"
			}
		}else{die "problem with triplets"}
		
	}else{
		#neither synonymous nor nonsynonymous
		$ptrAllData->{$chromosome}[$position]{syn}=-1;
	}
	return $ptrAllData;
}

sub translate{
	my ($nucleotide)=@_;
	my $t;
	if ($nucleotide eq "A"){$t = "T"}
	elsif($nucleotide eq "T"){$t = "A"}
	elsif($nucleotide eq "C"){$t = "G"}
	elsif($nucleotide eq "G"){$t = "C"}
	else{die "non translatable nucleotide letter $nucleotide."}
	return $t;
}

sub set_syn_nonsyn_nothing{
	my ($ptrAllData, $ptrGeneParts, $ptrCodons, $alphabet)=@_;
	
	foreach my $geneID (keys %{$ptrGeneParts}){
		foreach my $ptrPart (@{$ptrGeneParts->{$geneID}}){
			my $chromosome = $ptrPart->{chromosome};
			my $start = $ptrPart->{start};
			my $end = $ptrPart->{end};
			my $strand = $ptrPart->{strand};
			
			for (my $i=$start; $i<=$end; $i++){
				my $ptrPositionData = $ptrAllData->{$chromosome}[$i]; 
				$ptrAllData = set_syn_nonsyn_nothing_for_position($ptrAllData, $ptrPositionData, $chromosome, $i, $start, $end, $strand, $ptrCodons, $alphabet);	
	
			}
			
		}
	}

	
	return $ptrAllData; 
}


sub load_GPM_collect_genes_calculate_measures_and_overlaps{
	my (
		$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
		$POOL_SIZE,
		$SEGREGATING_SITES
	)=@_;
	
	my ($ptrAllData, $ptrGeneParts, $ptrWithdrawnInGtf) = load_all_data_GTF_PILEUP_MAUVEparsed(
											$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
											$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL
										);
	#$ptrGeneParts hash of arrays of well defined parts of the chromosom $ptrGeneParts->{geneID}[]{start, end, chromosome, strand}
#
#   decide, for each well defined position wheather the position is synonymous or nonsynonymous or neither of them.
#	
	my $ptrCodons = load_codon_table();
	
	my %hash; 
	foreach my $value (values %$ptrCodons){
		$hash{$value}=1
	}
	my $alphabet = join("", keys %hash);
	#print Dumper($alphabet);
	
	$ptrAllData = set_syn_nonsyn_nothing($ptrAllData, $ptrGeneParts, $ptrCodons, $alphabet); 
	my ($ptrGeneData, $ptrSplittedGenes) = store_genes_info($ptrGeneParts, $ptrAllData);	
	my $ptrGeneIDs = select_the_effectively_longest_transcript_for_each_gene($ptrGeneData);
	my $ptrOverlaps={}; my $ptrOverlappingGeneIDs={};
	($ptrOverlaps, $ptrOverlappingGeneIDs, $ptrGeneIDs) = find_overlaps_for_genes($ptrAllData, $ptrGeneIDs);

	$ptrGeneData = calculate_and_store_S_K_D_theta_avgCoverage($ptrGeneData, $ptrAllData, $ptrGeneIDs, $ptrGeneParts, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV, $SEGREGATING_SITES);
	$ptrGeneData = calculate_and_store_S_exp_and_var_for_genes($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData, $POOL_SIZE, $MIN_COUNT, $SEGREGATING_SITES);
		
	return ($ptrGeneIDs, $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf);	
}


#sub load_GPM_collect_genes_calculate_measures_and_overlaps_old_theta{
#	my (
#		$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
#		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,	
#		$POOL_SIZE,
#		$SEGREGATING_SITES
#	)=@_;	
#	
#	my ($ptrAllData, $ptrGeneParts, $ptrWithdrawnInGtf) = load_all_data_GTF_PILEUP_MAUVEparsed_old_theta(
#											$GTF_FILE, $WITHDRAWN_FILE, $PILEUP_FILE, $MAUVE_PARSED_FILE,
#											$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL
#										);
#	#$ptrGeneInfo hash of arrays of well defined parts of the chromosom $ptrGeneInfo->{geneID}[]{start, end, chromosome}
#	
#	my ($ptrGeneData, $ptrSplittedGenes) = store_genes_info($ptrGeneParts, $ptrAllData);	
#	my $ptrGeneIDs = select_the_effectively_longest_transcript_for_each_gene($ptrGeneData);
#	my $ptrOverlaps={}; my $ptrOverlappingGeneIDs={};
#	($ptrOverlaps, $ptrOverlappingGeneIDs, $ptrGeneIDs) = find_overlaps_for_genes($ptrAllData, $ptrGeneIDs);
#
#	calculate_and_store_measures_for_genes_old_theta($ptrGeneData, $ptrAllData, $ptrGeneIDs, $ptrGeneParts);
#	calculate_and_store_S_exp_and_var_for_genes($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData, $POOL_SIZE, $MIN_COUNT, $SEGREGATING_SITES);
#	
#	return ($ptrGeneIDs, $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf);	
#}



#############################################
####                                     ####
#### load previously calculated measures #### 
####                                     #### 
#############################################

# the slowest part of the script is loading data from GTF, PILEUP and MAUVE
# therefore, if you once have the measures calculated and stored the output somewhere, use it instead of loading GTF, PILEUP and MAUVE once more time.
# 
# loading files that are outputs of printing functions 

sub load_previously_calculated_measures_and_overlaps{
	my ($GENE_DATA_FILE, $OVERLAPPING_FILE)=@_;
	
	my $ptrGeneData = load_genes_measures($GENE_DATA_FILE);
	#print Dumper($ptrGeneData);
	my $ptrOverlappingGeneIDs = load_overlapping_genes($OVERLAPPING_FILE);
	
	return ($ptrGeneData, $ptrOverlappingGeneIDs);
}

sub load_genes_measures{
	#reading output from print_genes_measures()
	
	my ($FILE_NAME)=@_;
	my $ptrMeasures={};
	
	open fileHandle, "<", $FILE_NAME or die "Could not open file with measures for genes: $FILE_NAME.";
	while (my $line = <fileHandle>){
		next if ($line =~ m/^#/);
		chomp($line);
		my (
			$geneID, 
			$feature, 
			$chromosome, 
			$start, 
			$end, 
			$effectiveLengthOfGene, 
			$strand, 
			$avgCoverage, 
			$theta,   
			$Sfirst, 
			$SfirstExpectation, 
			$SfirstVariance, 
			$K, 
			$D
		) = split /\t/, $line;
		my $ptrTmpHash ={
			theta => $theta,
			K => $K,
			D => $D,
			Sfirst=> $Sfirst,
			SfirstExpectation=>$SfirstExpectation,
			SfirstVariance=>$SfirstVariance,
			feature=>$feature,
			chromosome => $chromosome,
			start=>$start,
			end=>$end,
			effectiveLengthOfGene => $effectiveLengthOfGene,
			strand=>$strand,
			avgCoverage => $avgCoverage, 			
		};
		
		$ptrMeasures->{$geneID}= $ptrTmpHash;
	}
	close fileHandle;
	return $ptrMeasures;
}

sub load_overlapping_genes{
#create a hash of overlapping geneIDs for each geneID
#reading output from print_overlapping_genes
	my ($FILE_NAME)=@_;
	
	my $ptrOverlappingGeneIDs={};
	
	open fileHandle, "<", $FILE_NAME or die "Could not open file with measures for genes: $FILE_NAME.";
	while (my $line = <fileHandle>){
		next if ($line =~ m/^#/);
		chomp($line);
		my ($geneIDkey, @geneIDval) = split /\t/, $line;
		$ptrOverlappingGeneIDs->{$geneIDkey} = \@geneIDval;
	}
	close fileHandle;
	return $ptrOverlappingGeneIDs;
}

sub load_reference_gene_list{
	my ($REFENRENCE_FILE)=@_;
	my $ptrReferenceGeneSet=[];
	
	open referenceHandle, "<", $REFENRENCE_FILE or die "Could not open file with a list of genes in reference set";
	while (my $line = <referenceHandle>){
		chomp($line);
		next if $line  =~ m/^#/;
		my @f = split "\t", $line; 
		push @{$ptrReferenceGeneSet}, $f[0];
	}
	close referenceHandle;
	return $ptrReferenceGeneSet;
}


###################################################################
####														   ####	
#### calculations of measures and finding overlaps among genes ####
####														   ####	
###################################################################

# subroutine create_actual_list_of_genes_on_position
#
# input: 
#	$ptrAllData->{chromosome}[position]={gtf, pileup, mauve info}
#	$ptrGeneIDs -- a list of genes, for each gene only the longest transcript is stored
# output:
#	 $ptrList->[] for a given position on chromosome output is a list of genes that overlaps on the position and are not transcripts of one gene 
#	
# used in: find_overlaps_for_genes
sub create_actual_list_of_genes_on_position{
	my ($ptrAllData, $chromosome, $position, $ptrGeneIDs)=@_;
	my $ptrList=[];
	
	foreach my $element ( @{$ptrAllData->{$chromosome}[$position]{geneID}} ){
		next unless contains($element, $ptrGeneIDs);
		push @$ptrList, $element;
	}
	return $ptrList;
}

# subroutine find_overlaps_for_genes
#
# input: 
#	$ptrAllData->{chromosome}[position]={gtf, pileup, mauve info}
# 	$ptrGeneIDs->[] -- a list of genes, for each gene only the longest transcript is stored
# output:
# 	$ptrOverlaps->{chromosome}[position]= [geneIDs that ovrlap on the position] 
# 	$ptrOverlappingGeneIDs->{geneID}= [geneIDs that overlaps somewhere with the hash key geneID]
#
# used in: load_GPM_collect_genes_calculate_measures_and_overlaps
sub find_overlaps_for_genes{
	my ($ptrAllData,$ptrGeneIDs)=@_;
	
	my $ptrOverlaps={};
	my $ptrOverlappingGeneIDs={}; 
	
	foreach my $chromosome (keys %$ptrAllData){
		for (my $i=0; $i< scalar @{$ptrAllData->{$chromosome}}; $i++){
			next unless defined($ptrAllData->{$chromosome}[$i]);
			my $ptrActListGeneIDs = create_actual_list_of_genes_on_position($ptrAllData, $chromosome, $i, $ptrGeneIDs);
			
			next unless (scalar @{$ptrActListGeneIDs} > 1);

			foreach my $geneID (@{$ptrActListGeneIDs}){
#update $ptrOverlaps
				push @{$ptrOverlaps->{$chromosome}[$i]}, $geneID;
#update $ptrOverlappingGeneIDs				
				my $ptrArray=$ptrOverlappingGeneIDs->{$geneID};
				foreach my $geneIDvalue (@{$ptrActListGeneIDs}){
					next if ($geneID eq $geneIDvalue);
					next if (contains($geneIDvalue, $ptrArray));
					push @{$ptrArray}, $geneIDvalue;						
				}	
				$ptrOverlappingGeneIDs->{$geneID}=$ptrArray;
			}
								
		}	
	}
	
	my $ptrUpdatedGeneIDsList=[];
	
	foreach my $geneID (@$ptrGeneIDs){
		next if exists $ptrOverlaps->{$geneID};
		push @{$ptrUpdatedGeneIDsList}, $geneID;
	}	
	
	return ($ptrOverlaps, $ptrOverlappingGeneIDs,$ptrUpdatedGeneIDsList);
}

# subroutine select_CDS_gene_well_defined_parts
# selects only the well defined (in mauve - not an unknown base - and pileup - covered) parts of gene
# output:
#	$ptrGeneIDs->{geneID} stores a ptrArray of info about continuous WELL DEFINED parts of a genome that are components of the gene
# input:
#	$ptrGeneIDsAllData->{geneID} stores a ptrArray of info about continuous parts of a genome that are components of the gene -- only CDS do not have to be covered enough etc.
#											ptrArray elements are {chromosome,start,end}
# used in: load_all_data_GTF_PILEUP_MAUVEparsed  
sub select_CDS_gene_well_defined_parts{
	my ($ptrAllData, $ptrGeneIDsAllData)=@_;
	
	my $ptrGeneParts={};
	
	foreach my $geneID (keys %$ptrGeneIDsAllData){
		foreach my $ptrGenePart (@{$ptrGeneIDsAllData->{$geneID}}){
			my $chromosome = $ptrGenePart->{chromosome};
			my $start =	$ptrGenePart->{start};
			my $end = $ptrGenePart->{end};
			my $strand = $ptrGenePart->{strand};

			my $i=$start;			
			my $partStart;
			my $partEnd;
			
			#iterate through the first undefined area
			while ( (!defined $ptrAllData->{$chromosome}[$i]) and ($i <= $end)){ $i+=1; }
			
			while ($i<= $end){ 
				$partStart=$i;
				
				#iterate through defined area
				while ( (defined $ptrAllData->{$chromosome}[$i]) and ($i<= $end) ){ $i+=1; }
			
				#shift one back from the first undefinef position
				$i-=1;
				$partEnd=$i;
			
				#write the part to $ptrGeneIDs
				if ($partStart<=$partEnd){
					push @{$ptrGeneParts->{$geneID}}, {chromosome=>$chromosome, start=>$partStart, end=>$partEnd, strand=>$strand};
				}
			
				#shift one position forward to the first undefined
				$i+=1;
				
				#iterate through undefined area
				if ($i<= $end){
					while ( (!defined $ptrAllData->{$chromosome}[$i]) and ($i <= $end)){ $i+=1; }
				}
			}
		}
	}

	return $ptrGeneParts;
}

# subroutine update_measure_and_length_sums
# updates iteratively many variables for gene measures calculation
# used in: calculate_measures_for_genes
sub update_measure_sums{
	my (
		$ptrAllData, 
		$chromosome, 
		$start, 
		$end, 
		$prevSumTheta, 
		$prevSumK, 
		$sumCoverage, 
		$prevSfirst, 
		$prevSsyn,
		$prevSnonsyn,
		$prevSothers, 
		$vec, 
		$MIN_COUNT, 
		$POOL_SIZE, 
		$MIN_COV, 
		$MAX_COV, 
		$SEGREGATING_SITES,
	) =@_;

	my $actSumTheta=0;
	my $actSumK=0;
	my $actSfirst=0;
	my $actSsyn=0;
	my $actSnonsyn=0;
	my $actSothers=0;
	
	for (my $i=$start; $i<=$end; $i++){
		my $ptrPositionData=$ptrAllData->{$chromosome}[$i];
		if ($ptrPositionData->{ispuresnp}){
			if($SEGREGATING_SITES =~ m/Corr/ ){
				$actSumTheta += $vec->calculate_measure("theta", [$ptrPositionData], 1);
			}else{
				$actSumTheta += 1/H_n_min_1_order_1($ptrPositionData->{eucov});
			}
			
			$actSfirst+=1;
			
			if ($ptrPositionData->{syn} == 1){	$actSsyn+=1	}
			elsif($ptrPositionData->{syn} == 0){	$actSnonsyn+=1	}
			elsif($ptrPositionData->{syn} == -1){	$actSothers+=1	}
			else{die "syn-nonsyn positions problem in update_measure_sums."}
			
			#print Dumper($ptrPositionData->{syn}, $actSfirst);		
		};
		
		# observed position is polymorphic
		if ($ptrPositionData->{firstRef} ne $ptrPositionData->{secondRef}){	$actSumK+=1	}
		
		$sumCoverage+= $ptrAllData->{$chromosome}[$i]{eucov};
	}
	
		
	my $updatedSumTheta= $prevSumTheta + $actSumTheta;
	my $updatedSumK=$prevSumK+$actSumK;
	my $updatedSfirst = $prevSfirst + $actSfirst;
	my $updatedSsyn = $prevSsyn + $actSsyn;
	my $updatedSnonsyn = $prevSnonsyn + $actSnonsyn;
	my $updatedSothers = $prevSothers + $actSothers;
	
	return (
		$updatedSumTheta, 
		$updatedSumK, 
		$sumCoverage, 
		$updatedSfirst,
		$updatedSsyn,
		$updatedSnonsyn,
		$updatedSothers
	);
}


sub update_min_max_bounds_of_gene{
	my ($min, $max, $start, $end)=@_;

	my $upMin = $min;
	my $upMax = $max;

	if($upMin==0){
		$upMin=$start;
	}else{
		if ($start<$upMin){
			$upMin=$start;
		}
	}

	if($upMax==0){
		$upMax=$end;
	}else{
		if($upMax<$end){
			$upMax=$end;
		}
	}
	return ($upMin,$upMax);
}

sub update_strand_values{
	my ($ptrStrands, $strand)=@_;
	
	if($strand eq "+"){
		$ptrStrands->{plus}+=1
	}elsif($strand eq "-"){
		$ptrStrands->{minus}+=1;
	};
	
}

sub set_strand_for_gene{
	my ($ptrStrands)=@_;
		
	my $strand;	 
	if (($ptrStrands->{plus} == 0)or($ptrStrands->{minus}==0)){
		if ($ptrStrands->{plus}==0){
			$strand="-"
		}else{
			$strand="+"
		}
	}else{
		$strand="+/-"
	};
	return $strand;		
}

sub update_list_of_chromosomes{
	my ($ptrChromosomes, $chromosome)=@_;
	
	$ptrChromosomes->{$chromosome} = 1;
}

sub update_length_of_gene{
	my ($length, $start, $end)=@_;
	
	my $upLength=$length;
	
	if ($start<=$end){
		$upLength+= $end - $start + 1;		
	}

	return $upLength;
}

# store_genes_info
# in: 
#	$ptrGeneParts->{geneID}[position]=={chromosome, start, end, strand}
#   $ptrAllData->{chromosome}[position]=={}
#chromosome,
#feature, 
#start, 
#end, 
#strand
sub store_genes_info{
	my ($ptrGeneParts, $ptrAllData)=@_;
	
	my $ptrGeneData={};
	my $ptrSplittedGenes={};
	
	
	my @geneIDlist;
	@geneIDlist = keys %$ptrGeneParts;
	
	foreach my $geneID (@geneIDlist){

		my $ptrChromosomes = {};
		my $ptrStrands = {plus=>0, minus=>0};
		
		my $min=0;
		my $max=0;
		my $length=0;
		
		foreach my $ptrPart (@{$ptrGeneParts->{$geneID}}){			
			my $chromosome = $ptrPart->{chromosome};
			my $start = $ptrPart->{start};
			my $end = $ptrPart->{end};
			my $strand = $ptrPart->{strand};
			$length=update_length_of_gene($length, $start, $end);
			update_list_of_chromosomes($ptrChromosomes, $chromosome);
			($min, $max) = update_min_max_bounds_of_gene($min, $max, $start, $end);
			update_strand_values($ptrStrands, $strand);
		};
			
		if (scalar keys %$ptrChromosomes == 1){
		# gene is situated only on one chromosome 	
			my $chromosome = (keys %$ptrChromosomes)[0];
			my $strand = set_strand_for_gene($ptrStrands);
			my $feature;
			if ("C" eq $ptrAllData->{$chromosome}[$min]{feature}){$feature="CDS"};
			$ptrGeneData->{$geneID}={
				chromosome => $chromosome,
				start=>$min,
				end=>$max,
				strand=>$strand,
				feature=>$feature,
				effectiveLengthOfGene=>$length
			}
		}else{
 		# gene is splitted 	
			my @chromosomesArray = keys %$ptrChromosomes;
			$ptrSplittedGenes->{$geneID}=\@chromosomesArray;	
		}
	}
	return ($ptrGeneData, $ptrSplittedGenes);
}


# subroutine calculate_and_store_S_K_D_theta_avgCoverage
# calculates the measures we need for HKA test and stores some additional info for further analysis
# output:
# 	$ptrMeasures: $ptrMeasures->{geneID} stores a hash of info about gene
#   $ptrMeasuresSplittedGenes: the same info as in $ptrMeasure for genes that are splitted across more than one gene     
# info: 
#		feature -- so far only CDS 
#		chromosome
#		first chr position
#		last chr position
#		strand		
#		theta, 
#		theta correction factor,
#		theta correction sqr factor,
#		S
#		exp[S]
#		var[S]
#		K

sub calculate_and_store_S_K_D_theta_avgCoverage{
	my ($ptrGeneData, $ptrAllData, $ptrGeneIDs, $ptrGeneParts, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV, $SEGREGATING_SITES) = @_;
	
	my $sumTheta;
	my $sumK;
	my $length;
	my $sumCoverage;
	my $Sfirst;
	my $Ssyn;
	my $Snonsyn;
	my $Sothers;
		
	my $vec;
		
	#if($SEGREGATING_SITES =~ m/Cn/){
	#old way of theta calculation, per site sum of 1/Hn-1	
	#}elsif($SEGREGATING_SITES =~ m/Corr/){
	#corrected theta calculation, per site correction
	
	if($SEGREGATING_SITES =~ m/Corr/){
		$vec = VarianceExactCorrection->new($POOL_SIZE, $MIN_COUNT);
	}	
	 	
	foreach my $geneID (@$ptrGeneIDs){
		$sumTheta=0;
		$sumK=0;
		$length= $ptrGeneData->{$geneID}{effectiveLengthOfGene};
		$sumCoverage=0;
		$Sfirst=0;
		$Ssyn=0;
		$Snonsyn=0;
		$Sothers=0;
		
		foreach my $ptrPart (@{$ptrGeneParts->{$geneID}}){			
			my $chromosome = $ptrPart->{chromosome};
			my $start = $ptrPart->{start};
			my $end = $ptrPart->{end}; 
				
			($sumTheta, $sumK, $sumCoverage, $Sfirst, $Ssyn, $Snonsyn, $Sothers) = update_measure_sums(	$ptrAllData, 
																										$chromosome, 
																										$start, 
																										$end, 
																										$sumTheta,  
																										$sumK, 
																										$sumCoverage, 
																										$Sfirst,
																										$Ssyn,
																										$Snonsyn,
																										$Sothers, 
																										$vec, 
																										$MIN_COUNT, 
																										$POOL_SIZE, 
																										$MIN_COV, 
																										$MAX_COV,
																										$SEGREGATING_SITES
																					);
		};
		
		my $theta=0;
		my $K=0;
		my $avgCoverage=0;
		
		if ($length != 0){ 
			$avgCoverage = $sumCoverage / $length;
			my $avgCoverageRound = int($avgCoverage+0.5);
			$theta = $sumTheta / $length;
			$K = $sumK / $length;	
		}	
	
		$ptrGeneData->{$geneID}{theta}=$theta;
		$ptrGeneData->{$geneID}{K}=$K;
		$ptrGeneData->{$geneID}{D}=$sumK;
		$ptrGeneData->{$geneID}{Sfirst}=$Sfirst;
		$ptrGeneData->{$geneID}{Ssyn}=$Ssyn;
		$ptrGeneData->{$geneID}{Snonsyn}=$Snonsyn;
		$ptrGeneData->{$geneID}{Sothers}=$Sothers;
		$ptrGeneData->{$geneID}{avgCoverage}=$avgCoverage;
		
		#print Dumper($Sfirst, $Ssyn, $Snonsyn, $Sothers);
	};	
	return $ptrGeneData;	
};









sub calculate_and_store_S_exp_and_var_for_genes{
	my ($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData, $POOL_SIZE, $MIN_COUNT, $SEGREGATING_SITES)=@_;

	if ($SEGREGATING_SITES eq "avgCn"){
		$ptrGeneData = calculate_and_store_S_exp_and_var_for_genes_avg_H_n($ptrGeneData, $ptrGeneIDs);
	}elsif($SEGREGATING_SITES eq "perSiteCn"){
		$ptrGeneData = calculate_and_store_S_exp_and_var_for_genes_H_n_per_site($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData);
	}elsif($SEGREGATING_SITES eq "avgCorr"){
		$ptrGeneData = calculate_and_store_S_exp_and_var_for_genes_avg_corr($ptrGeneData, $ptrGeneIDs, $POOL_SIZE, $MIN_COUNT);
	}elsif($SEGREGATING_SITES eq "perSiteCorr"){
		$ptrGeneData = calculate_and_store_S_exp_and_var_for_genes_corr_per_site($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData, $POOL_SIZE, $MIN_COUNT)
	};
	
	return $ptrGeneData;
	
}

sub calculate_and_store_S_exp_and_var_for_genes_corr_per_site{
	my ($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData, $POOL_SIZE, $MIN_COUNT)=@_;
	
	foreach my $geneID (@$ptrGeneIDs){
		my $theta = $ptrGeneData->{$geneID}{theta};
		my $length= $ptrGeneData->{$geneID}{effectiveLengthOfGene};
		my $Sfirst = $ptrGeneData->{$geneID}{Sfirst};
		
		my $SfirstExpectation=0;
		my $SfirstVariance=0;
		
		my $div_buffer = get_thetadiv_buffer();
		my $div_buffer_sqr = get_thetadiv_buffer_sqr();
		

		foreach my $ptrGenePart (@{$ptrGeneParts->{$geneID}}){
			my $chromosome = $ptrGenePart->{chromosome};
			my $start = $ptrGenePart->{start};
			my $end = $ptrGenePart->{end}; 
			
			#update SfirstExpectation according to positions in $ptrGenePart 
			for (my $i=$start; $i<=$end; $i++){
				my $coverage = $ptrAllData->{$chromosome}[$i]{eucov};
		
				my $correction = $div_buffer->($MIN_COUNT, $POOL_SIZE, $coverage);
				#my $correction_sqr = $div_buffer_sqr->($MIN_COUNT, $POOL_SIZE, $coverage);		
				my $correction_sqr = H_n_min_1_order_2($coverage);
				
				$SfirstExpectation+=$theta*$correction;
				$SfirstVariance+=$correction_sqr * $length * ($theta)**2;
			}
		}
		$SfirstVariance += $SfirstExpectation;
		
		$ptrGeneData->{$geneID}{SfirstExpectation}=$SfirstExpectation;
		$ptrGeneData->{$geneID}{SfirstVariance}=$SfirstVariance;				
	}

	return $ptrGeneData;
	
}

sub calculate_and_store_S_exp_and_var_for_genes_avg_corr{
	my ($ptrGeneData, $ptrGeneIDs, $POOL_SIZE, $MIN_COUNT)=@_;

	foreach my $geneID (@$ptrGeneIDs){
		my $theta = $ptrGeneData->{$geneID}{theta};
		my $length= $ptrGeneData->{$geneID}{effectiveLengthOfGene};
		
		my $avgCoverage =  $ptrGeneData->{$geneID}{avgCoverage};
		my $roundAvgCoverage = int($avgCoverage+0.5);

		my $div_buffer = get_thetadiv_buffer();
		my $div_buffer_sqr = get_thetadiv_buffer_sqr();
		
		my $correction = $div_buffer->($MIN_COUNT, $POOL_SIZE, $roundAvgCoverage);
		#my $correction_sqr = $div_buffer_sqr->($MIN_COUNT, $POOL_SIZE, $roundAvgCoverage);
		my $correction_sqr = H_n_min_1_order_2($roundAvgCoverage);
		
		my $SfirstExpectation = $theta*$correction*$length;
		my $SfirstVariance = $SfirstExpectation + $correction_sqr * ($theta * $length)**2;	 
		
		$ptrGeneData->{$geneID}{SfirstExpectation}=$SfirstExpectation;
		$ptrGeneData->{$geneID}{SfirstVariance}=$SfirstVariance;				
	}

	return $ptrGeneData;
	
}


sub calculate_and_store_S_exp_and_var_for_genes_avg_H_n{
	my ($ptrGeneData, $ptrGeneIDs)=@_;
	
	foreach my $geneID (@$ptrGeneIDs){
		my $theta = $ptrGeneData->{$geneID}{theta};
		my $length= $ptrGeneData->{$geneID}{effectiveLengthOfGene};
		my $Sfirst = $ptrGeneData->{$geneID}{Sfirst};
		my $avgCoverage =  $ptrGeneData->{$geneID}{avgCoverage};
		my $roundAvgCoverage = int($avgCoverage+0.5);
		
		my $SfirstExpectation=$theta * $length * H_n_min_1_order_1($roundAvgCoverage);
		my $SfirstVariance=$SfirstExpectation + H_n_min_1_order_2($roundAvgCoverage) * ($theta * $length)**2;
		
		$ptrGeneData->{$geneID}{SfirstExpectation}=$SfirstExpectation;
		$ptrGeneData->{$geneID}{SfirstVariance}=$SfirstVariance;				
	}
	
	return $ptrGeneData;
}

sub calculate_and_store_S_exp_and_var_for_genes_H_n_per_site{
	my ($ptrGeneData, $ptrGeneIDs, $ptrGeneParts, $ptrAllData)=@_;
	
	foreach my $geneID (@$ptrGeneIDs){
		my $theta = $ptrGeneData->{$geneID}{theta};
		my $length= $ptrGeneData->{$geneID}{effectiveLengthOfGene};
		my $Sfirst = $ptrGeneData->{$geneID}{Sfirst};
		
		my $SfirstExpectation=0;
		my $SfirstVariance=0;
		
		foreach my $ptrGenePart (@{$ptrGeneParts->{$geneID}}){
			my $chromosome = $ptrGenePart->{chromosome};
			my $start = $ptrGenePart->{start};
			my $end = $ptrGenePart->{end}; 
			
			#update SfirstExpectation according to positions in $ptrGenePart 
			for (my $i=$start; $i<=$end; $i++){
				my $coverage = $ptrAllData->{$chromosome}[$i]{eucov};
				$SfirstExpectation+=$theta*H_n_min_1_order_1($coverage);
				$SfirstVariance+=H_n_min_1_order_2($coverage) * $length * ($theta)**2;
			}
		}
		
		$SfirstVariance += $SfirstExpectation;
		$ptrGeneData->{$geneID}{SfirstExpectation}=$SfirstExpectation;
		$ptrGeneData->{$geneID}{SfirstVariance}=$SfirstVariance;								
	}	
	
	return $ptrGeneData;
}

#sub calculate_and_store_measures_for_genes_old_theta{
#	my ($ptrGeneData, $ptrAllData, $ptrGeneIDs, $ptrGeneParts) = @_;
#
#	foreach my $geneID (@$ptrGeneIDs){
#		my $sumTheta=0;
#		my $sumK=0;
#		my $length= $ptrGeneData->{$geneID}{effectiveLengthOfGene};
#		my $sumCoverage=0;
#		my $Sfirst=0;
#		
#		foreach my $ptrPart (@{$ptrGeneParts->{$geneID}}){			
#			my $chromosome = $ptrPart->{chromosome};
#			my $start = $ptrPart->{start};
#			my $end = $ptrPart->{end}; 
#			
#			($sumTheta, $sumK, $sumCoverage, $Sfirst) = update_measure_sums_old_theta($ptrAllData, 
#																			$chromosome, 
#																			$start, 
#																			$end, 
#																			$sumTheta,  
#																			$sumK, 
#																			$sumCoverage, 
#																			$Sfirst, 
#																		);
#		};
#		
#		my $theta=0;
#		my $K=0;
#		my $avgCoverage =0;
#		
#		if ($length != 0){ 
#			$avgCoverage = $sumCoverage / $length;
#			my $avgCoverageRound = int($avgCoverage+0.5);
#			$theta = $sumTheta / $length;
#			$K = $sumK / $length;	
#		}	
#	
#		$ptrGeneData->{$geneID}{theta}=$theta;
#		$ptrGeneData->{$geneID}{K}=$K;
#		$ptrGeneData->{$geneID}{D}=$sumK;
#		$ptrGeneData->{$geneID}{Sfirst}=$Sfirst;
#		$ptrGeneData->{$geneID}{avgCoverage}=$avgCoverage;
#	};
#};





sub H_n_min_1_order_1{
	my ($n)=@_;
	my $sum=0;
	
	for (my $i=1; $i<$n; $i++){
		$sum+=1/$i;
	}
	return $sum;
}

sub H_n_min_1_order_2{
	my ($n)=@_;
	
	my $sum=0;
	
	for (my $i=1; $i<$n; $i++){
		$sum+=1/$i**2;
	}
	return $sum;	
}
# select_the_effectively_longest_transcript_for_each_gene
# 
# in:
#	$ptrGeneData->{geneID}{chromosome, start, end, strand, feature, effectiveLengthOfGene}
# out: 
#	$ptrGeneIDs->[position]==geneID
# used in: 
#	load_GPM_collect_genes_calculate_measures_and_overlaps
sub select_the_effectively_longest_transcript_for_each_gene{
	my ($ptrGeneData)=@_;
	
	my $ptrGeneGroups={};
	
	my @geneIDlist = keys %$ptrGeneData;
	
	foreach my $geneID (@geneIDlist){
		my ($geneGr, $transcript) = split /-/,$geneID;
		
		my $newEffLength = $ptrGeneData->{$geneID}{effectiveLengthOfGene};
		
		if (exists($ptrGeneGroups->{$geneGr})){
			if ($ptrGeneGroups->{$geneGr}{effectiveLengthOfGene} < $newEffLength){
				$ptrGeneGroups->{$geneGr}{geneID} = $geneID;
				$ptrGeneGroups->{$geneGr}{effectiveLengthOfGene}=$newEffLength;
			}
		}else{
			$ptrGeneGroups->{$geneGr}{geneID}= $geneID;
			$ptrGeneGroups->{$geneGr}{effectiveLengthOfGene}=$newEffLength;
		} 
	}

	my $ptrGeneIDs=[];
	foreach my $geneID (@geneIDlist){
		my ($geneGr, $transcript) = split /-/, $geneID;
		next unless ($ptrGeneGroups->{$geneGr}{geneID} eq $geneID);
		push @$ptrGeneIDs, $geneID;
	}

	return $ptrGeneIDs;
}
		
############################################################       
####                                                    #### 
#### calculating HKAtest statistic and needed variables ####
####                                                    ####
############################################################

#
# random selection of $COUNT genes
#

sub handling_with_overlapping_and_short_and_zero_variance_genes{
	#returns measures for selected non-overrlapping and long-enough genes
	my ($ptrGeneData, $ptrOverlappingGeneIDs, $MIN_LENGTH, $ptrGeneIDs)=@_;
	
	my $ptrUpdatedGeneIDs = [];
	my $ptrRemovedGenes = {};
	
	($ptrUpdatedGeneIDs, $ptrRemovedGenes) = create_list_of_nonoverlapping_genes($ptrGeneData, $ptrOverlappingGeneIDs, $ptrGeneIDs);
	($ptrUpdatedGeneIDs, $ptrRemovedGenes) = create_list_of_long_enough_genes($ptrGeneData, $ptrUpdatedGeneIDs, $ptrRemovedGenes, $MIN_LENGTH);
	($ptrUpdatedGeneIDs, $ptrRemovedGenes) = create_list_of_genes_with_non_zero_S_variance($ptrGeneData, $ptrUpdatedGeneIDs, $ptrRemovedGenes);
	
	return ($ptrUpdatedGeneIDs, $ptrRemovedGenes);
}

sub create_list_of_nonoverlapping_genes{
	my ($ptrGeneData, $ptrOverlappingGeneIDs, $ptrGeneIDs)=@_;
	
	my $ptrUpdatedGeneIDs=[];
	my $ptrRemovedGenes ={};
	
	my @listOfAllGenes;
	
	if (defined($ptrGeneIDs)){
		@listOfAllGenes = @$ptrGeneIDs;
	} else {
		@listOfAllGenes = keys %$ptrGeneData;
	}
	
	foreach my $geneID (sort @listOfAllGenes){
		if (defined($ptrOverlappingGeneIDs->{$geneID})){
			$ptrRemovedGenes->{$geneID}= "overlaps with genes: @{$ptrOverlappingGeneIDs->{$geneID}}.";
		}else{
			push @$ptrUpdatedGeneIDs, $geneID;
		}
	}	
	return ($ptrUpdatedGeneIDs, $ptrRemovedGenes);
};

sub create_list_of_long_enough_genes{
	my ($ptrGeneData, $ptrGeneIDs, $ptrRemovedGenes, $MIN_LENGTH)=@_;
	
	my $ptrUpdatedGeneIDs=[];
	my $ptrRemovedUpdate = $ptrRemovedGenes;
	
	foreach my $geneID (@{$ptrGeneIDs}){
		#print Dumper($ptrGeneData->{$geneID}{effectiveLengthOfGene});
		if ($ptrGeneData->{$geneID}{effectiveLengthOfGene} < $MIN_LENGTH){
			$ptrRemovedUpdate->{$geneID}="too short gene, length: ".$ptrGeneData->{$geneID}{effectiveLengthOfGene}.".";
		}else{
			push @$ptrUpdatedGeneIDs, $geneID;
		}
	}
	return ($ptrUpdatedGeneIDs, $ptrRemovedUpdate);
};

sub create_list_of_genes_with_non_zero_S_variance{
	my ($ptrGeneData, $ptrGeneIDs, $ptrRemovedGenes)=@_;
	
	open errHandle, ">",  "problem.log" or die "Could not open log file.";
	my $ptrUpdatedGeneIDs=[];
	my $ptrRemovedUpdate=$ptrRemovedGenes;
	#print Dumper($ptrGeneData);
	foreach my $geneID (@{$ptrGeneIDs}){
		#print $geneID, "\n";
		if (!defined($ptrGeneData->{$geneID}{Sfirst})){
			
			print errHandle $geneID, "\n"; 
		}else{
		
		if ($ptrGeneData->{$geneID}{Sfirst}==0){
			$ptrRemovedUpdate->{$geneID}="gene with no segregating sites.";
		}else{
			push @$ptrUpdatedGeneIDs, $geneID;
		}
		
		}
	}
	close errHandle;
	return ($ptrUpdatedGeneIDs, $ptrRemovedUpdate);
}

sub randomly_choose_genes{
	my ($ptrMeasures, $ptrGeneIDs, $COUNT) = @_;

	my $ptrRandomlySelectedGenes = [];
	
	my $numberOfGenes = scalar @$ptrGeneIDs;
	my @geneList = sort @$ptrGeneIDs;
	
	if ($numberOfGenes < $COUNT){
		die "random selection of genes failed: not so many genes as $COUNT"
	}else{
		srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);	
		#srand(10);	
	
		my $pos=0;
		my $remains = $COUNT;
		my $remainsAll = $numberOfGenes;
		my $probability = $COUNT / $numberOfGenes;
	
		while ( ($remains>0) and ($probability>0) and ($pos<$numberOfGenes)){
			my $r = rand();
			#print Dumper($pos,$probability,$r, $remains, $remainsAll);
			if ($r<=$probability){
				push @$ptrRandomlySelectedGenes, $geneList[$pos];
				$remains-=1;
			};
			$pos+=1;
			$remainsAll-=1;
			if ($remainsAll!=0){$probability = $remains / $remainsAll;}
		}	
		
	}	
	return $ptrRandomlySelectedGenes;
}


#
# calculate the remains measures, that depends on a sample group
#  T, E(D), Var(D)
#

sub calculate_and_store_T{
	my ($ptrGeneData, $ptrRandomGenesList, $selectedGeneID)=@_;
	
	#calculate
	my $sumD = 0;
	my $sumLengthTheta= 0;
	
	if (defined($selectedGeneID)){
		$sumD = $ptrGeneData->{$selectedGeneID}{D};
		$sumLengthTheta=($ptrGeneData->{$selectedGeneID}{effectiveLengthOfGene}) * ($ptrGeneData->{$selectedGeneID}{theta});
	}
	
	foreach my $geneID (@$ptrRandomGenesList){
		$sumD += $ptrGeneData->{$geneID}{D};
		$sumLengthTheta += ($ptrGeneData->{$geneID}{effectiveLengthOfGene}) * ($ptrGeneData->{$geneID}{theta}); 
	}
	
	my $T = $sumD/$sumLengthTheta - 1;
	
	if (defined($selectedGeneID)){
		$ptrGeneData->{$selectedGeneID}{T}=$T;
	}else{
		foreach my $geneID (@$ptrRandomGenesList){
			$ptrGeneData->{$geneID}{T}=$T;
		}
	}
	
	if (defined($selectedGeneID)){
		foreach my $geneID (@{$ptrRandomGenesList}){
			$ptrGeneData->{$geneID}{$selectedGeneID}{T}=$T;
		}
	}
}

sub calculate_and_store_D_expectations_and_variance{
	my ($ptrGeneData, $ptrRandomGenesList, $selectedGeneID)=@_;
	
	foreach my $geneID (@$ptrRandomGenesList){
		my $Dexpectation = $ptrGeneData->{$geneID}{theta} * $ptrGeneData->{$geneID}{effectiveLengthOfGene} * ($ptrGeneData->{$geneID}{$selectedGeneID}{T} + 1) ;
		
		$ptrGeneData->{$geneID}{$selectedGeneID}{Dexpectation}=$Dexpectation;
		$ptrGeneData->{$geneID}{$selectedGeneID}{Dvariance}= $Dexpectation + ($ptrGeneData->{$geneID}{theta}*$ptrGeneData->{$geneID}{effectiveLengthOfGene})**2;
	}
	
	my $Dexp = $ptrGeneData->{$selectedGeneID}{theta} * ($ptrGeneData->{$selectedGeneID}{T} + 1) * $ptrGeneData->{$selectedGeneID}{effectiveLengthOfGene};
		
	$ptrGeneData->{$selectedGeneID}{Dexpectation}=$Dexp;
	$ptrGeneData->{$selectedGeneID}{Dvariance}= $Dexp + ($ptrGeneData->{$selectedGeneID}{theta}*$ptrGeneData->{$selectedGeneID}{effectiveLengthOfGene})**2;
}

sub calculate_and_store_chi_square_statistic_for_one_gene{
#HKA chi square statistic from paper HKA 1987, without S_i^B sum
# the chi_square_stat is calculated in the following way:
#	
#	chi_square_stat = sum_{i=1}^L  { (S_i^A - E[S_i^A])^2/Var[S_i^A] + (D_i - E[D_i])^2/Var[D_i] }
#
	my ($ptrGeneData, $ptrReferenceSetOfGenes, $selectedGeneID, $degreesOfFreedom)=@_;
	
	#my $nullVarianceD;
	
	#$nullVarianceD = ($ptrGeneData->{$selectedGeneID}{Dvariance} == 0);
	
	my $sum=0;
	# CHANGED
	$sum+= ($ptrGeneData->{$selectedGeneID}{Sfirst} - $ptrGeneData->{$selectedGeneID}{SfirstExpectation})**2 / $ptrGeneData->{$selectedGeneID}{SfirstVariance};

	#if (!$nullVarianceD){
		$sum+= 	($ptrGeneData->{$selectedGeneID}{D} - $ptrGeneData->{$selectedGeneID}{Dexpectation})**2/ $ptrGeneData->{$selectedGeneID}{Dvariance};
	#}
	
	foreach my $geneID (@$ptrReferenceSetOfGenes){

		#$nullVarianceD = ($ptrGeneData->{$geneID}{$selectedGeneID}{Dvariance} == 0);
		#CHANGED
		$sum+= ($ptrGeneData->{$geneID}{Sfirst} - $ptrGeneData->{$geneID}{SfirstExpectation})**2 / $ptrGeneData->{$geneID}{SfirstVariance};
	
		#if (!$nullVarianceD){
			$sum+= ($ptrGeneData->{$geneID}{D} - $ptrGeneData->{$geneID}{$selectedGeneID}{Dexpectation})**2/ $ptrGeneData->{$geneID}{$selectedGeneID}{Dvariance};
		#}
	}
	
	$ptrGeneData->{$selectedGeneID}{chiSquare}=$sum;
	$ptrGeneData->{$selectedGeneID}{degreesOfFreedom}=$degreesOfFreedom;
	$ptrGeneData->{$selectedGeneID}{pValue}=chisqrprob($degreesOfFreedom, $sum);
}

sub calculate_and_store_chi_square_statistic_for_each_gene{
	my ($ptrGeneData, $ptrReferenceSetOfGenes, $ptrGeneIDs)=@_;
	
	calculate_and_store_T($ptrGeneData, $ptrReferenceSetOfGenes);
	
	foreach my $selectedGeneID (@{$ptrGeneIDs}){
		next if contains($selectedGeneID, $ptrReferenceSetOfGenes);
		
		#calculate
		calculate_and_store_T($ptrGeneData, $ptrReferenceSetOfGenes, $selectedGeneID);
		calculate_and_store_D_expectations_and_variance($ptrGeneData, $ptrReferenceSetOfGenes, $selectedGeneID);
		my $degreesOfFreedom = scalar @{$ptrReferenceSetOfGenes}; 
		calculate_and_store_chi_square_statistic_for_one_gene($ptrGeneData, $ptrReferenceSetOfGenes, $selectedGeneID, $degreesOfFreedom);
		 
	}
	
}
 
# load genesData and overlapping genes
# handling_with_overlapping_and_short_genes
# calculate_and_store_chi_square_for_each_gene

# print info about all other genes with chi square stat
# print info about randomly selected genes

# TODO: print params, timestamps

#TODO: add commented headers to print functions with runned file names, parameters and columns descriptors


 
 
 
##################################### 
##                                 ##  
## printing functions of all kinds ##
##                                 ## 
##################################### 
 
sub print_first_round_of_data{
	my ($ptrGeneIDs, $ptrGeneData, $ptrSplittedGenes, $ptrOverlaps, $ptrOverlappingGeneIDs, $ptrWithdrawnInGtf, $filehandleData, $filehandleOverlappingGenes, $filehandleLog)=@_;
	
	print_selected_genes_without_T($ptrGeneIDs, $ptrGeneData, $filehandleData);
	print_overlapping_genes($ptrOverlappingGeneIDs, $filehandleOverlappingGenes);
	print_log_multiple_geneIDs_positions($ptrOverlaps, $filehandleLog);
	print_log_splitted_genes($ptrSplittedGenes, $filehandleLog);
	print_log_withdrawn($ptrWithdrawnInGtf, $filehandleLog);
}

sub print_log_withdrawn{
	my ($ptrWithdrawnInGtf, $filehandle)=@_;
	
	print $filehandle "withdrawn genes that occured in input gtf file:\n";
	foreach my $geneID (keys %$ptrWithdrawnInGtf){
		print $filehandle $geneID, "\n";
	}
}

sub print_log_removed_genes{
	my ($ptrRemovedGenes, $filehandle)=@_;
	
	print $filehandle "genes removed from further analysis and reasons why they were removed:\n";
	foreach my $geneID (keys %$ptrRemovedGenes){
		print $filehandle $geneID."\t".$ptrRemovedGenes->{$geneID}."\n";
	}
}

sub array_to_tab_separ_string{
	my ($ptrArray)=@_;
	
	my $str = "";
	if (scalar @$ptrArray == 0){
		return $str;
	}else{
		$str=$ptrArray->[0];
		for (my $i=1; $i< scalar @{$ptrArray}; $i++){
			$str=$str."\t".$ptrArray->[$i];
		}
		return $str;
	}
}

sub print_overlapping_genes{
	my ($ptrOverlappingGeneIDs, $filehandle)=@_;
	
	if (scalar(keys %{$ptrOverlappingGeneIDs}) > 0){
		print $filehandle "#The following different genes are overlapping each other.\n";
		#TODO: fix this problem,... select a nonoverlapping subset or so,...
		print $filehandle "#So far the genes will be excluded from a further analysis. \n";
		print $filehandle "#geneID\tarray_of_genes_that_overlap_with_the_geneID\n";
		foreach my $geneIDkey (keys %{$ptrOverlappingGeneIDs}){
			print $filehandle "$geneIDkey\t".array_to_tab_separ_string($ptrOverlappingGeneIDs->{$geneIDkey})."\n";
		}
	}else{
		print $filehandle "#There are no overlapping genes.\n";
	}
	
} 

sub print_log_multiple_geneIDs_positions{
	my ($ptrOverlaps, $filehandle)=@_;
		
	if (scalar(keys %{$ptrOverlaps})>0){
		print $filehandle "Positions on which different genes overlaps:\n";
		foreach my $chromosome (keys %$ptrOverlaps){
			for (my $i=0; $i< scalar @{$ptrOverlaps->{$chromosome}}; $i++){
				next unless defined($ptrOverlaps->{$chromosome}[$i]);
				next unless (scalar @{$ptrOverlaps->{$chromosome}[$i]} > 1);
				print $filehandle $chromosome, "\t", $i, "\t", "@{$ptrOverlaps->{$chromosome}[$i]}", "\n";
		}
	}
	}else{
		print $filehandle "There are no gene overlapping places in the the given data.";
	}
}

sub print_log_splitted_genes{
	my ($ptrSplittedGenes, $filehandle)=@_;

	if ( scalar(keys %$ptrSplittedGenes) > 0){
		print $filehandle "The following genes are splitted across more then one contig/chromosome and will be excluded from a further analysis:\n";
		foreach my $geneID (keys %$ptrSplittedGenes){
			print $filehandle "#geneID\t array of chromosomes\n";
			print $filehandle $geneID,"\t",(keys %{$ptrSplittedGenes->{$geneID}{chromosomes}}), "\n";
		}
	}else{
		print $filehandle "There are no genes splitted across more then one contig/chromosome.\n";
	}
}

sub print_selected_genes_without_T{
	my ($ptrSelectedGenes, $ptrGeneData, $filehandle)=@_;
	
#	print $filehandle "# S -- number of segregating sites within one spieces in the gene\n";
#	print $filehandle "# D -- number of divergent sites between spieces in the gene\n";
	my $header="";
	my $line="";
	
	$header="#geneID"."\t".
			"feature"."\t".
			"chromosome"."\t".
			"start_position"."\t".
			"end_position"."\t".
			"effective_length_of_gene"."\t".
			"strand"."\t".
			"average_coverage"."\t".
			"theta"."\t".
			"S_first"."\t".
			"S_first_syn"."\t".
			"S_first_nonsyn"."\t".
			"S_first_others"."\t".
			"S_first_expectation"."\t".
			"S_first_variance"."\t".
			"K"."\t".
			"D"."\n";
		
	print $filehandle $header;
	
	foreach my $selectedGeneID (@{$ptrSelectedGenes}){
		$line= 
			$selectedGeneID."\t". 
			$ptrGeneData->{$selectedGeneID}{feature}."\t".
			$ptrGeneData->{$selectedGeneID}{chromosome}."\t".
			$ptrGeneData->{$selectedGeneID}{start}."\t".
			$ptrGeneData->{$selectedGeneID}{end}."\t".
			$ptrGeneData->{$selectedGeneID}{effectiveLengthOfGene}."\t".
			$ptrGeneData->{$selectedGeneID}{strand}."\t".
			$ptrGeneData->{$selectedGeneID}{avgCoverage}."\t".
			$ptrGeneData->{$selectedGeneID}{theta}."\t".
			$ptrGeneData->{$selectedGeneID}{Sfirst}."\t".
			$ptrGeneData->{$selectedGeneID}{Ssyn}."\t".
			$ptrGeneData->{$selectedGeneID}{Snonsyn}."\t".
			$ptrGeneData->{$selectedGeneID}{Sothers}."\t".
			$ptrGeneData->{$selectedGeneID}{SfirstExpectation}."\t".
			$ptrGeneData->{$selectedGeneID}{SfirstVariance}."\t".
			$ptrGeneData->{$selectedGeneID}{K}."\t".
			$ptrGeneData->{$selectedGeneID}{D};
		
		$line=$line."\n";
		print $filehandle $line;			
	}
}

sub print_selected_genes_with_all_statistics{
	my ($ptrSelectedGenes, $ptrGeneData, $filehandle)=@_;
	
#	print $filehandle "# S -- number of segregating sites within one spieces in the gene\n";
#	print $filehandle "# D -- number of divergent sites between spieces in the gene\n";
	my $header="";
	my $line="";
	
	$header="#geneID"."\t".
			"chi_square"."\t".
			"degrees_of_freedom"."\t".
			"p-value"."\t".
			
			"average_coverage"."\t".
			"effective_length_of_gene"."\t".
			"chromosome"."\t".
			"start_position"."\t".
			"end_position"."\t".
			"feature"."\t".
			"strand"."\t".
			
			"theta"."\t".
			"S_first"."\t".			
			"S_first_expectation"."\t".
			"T"."\t".
			"D"."\t".
			"D_expectation"."\t".
			"K"."\t".
			"S_first_variance"."\t".
			"D_variance"."\n";
	
	print $filehandle $header;
	
	foreach my $selectedGeneID (@{$ptrSelectedGenes}){
		$line= 
			$selectedGeneID."\t". 
			$ptrGeneData->{$selectedGeneID}{chiSquare}."\t".
			$ptrGeneData->{$selectedGeneID}{degreesOfFreedom}."\t".
			$ptrGeneData->{$selectedGeneID}{pValue}."\t".
			
			$ptrGeneData->{$selectedGeneID}{avgCoverage}."\t".
			$ptrGeneData->{$selectedGeneID}{effectiveLengthOfGene}."\t".
			$ptrGeneData->{$selectedGeneID}{chromosome}."\t".
			$ptrGeneData->{$selectedGeneID}{start}."\t".
			$ptrGeneData->{$selectedGeneID}{end}."\t".
			$ptrGeneData->{$selectedGeneID}{feature}."\t".
			$ptrGeneData->{$selectedGeneID}{strand}."\t".
			
			$ptrGeneData->{$selectedGeneID}{theta}."\t".
			$ptrGeneData->{$selectedGeneID}{Sfirst}."\t".
			$ptrGeneData->{$selectedGeneID}{SfirstExpectation}."\t".
			$ptrGeneData->{$selectedGeneID}{T}."\t".
			$ptrGeneData->{$selectedGeneID}{D}."\t".
			$ptrGeneData->{$selectedGeneID}{Dexpectation}."\t".
			$ptrGeneData->{$selectedGeneID}{K}."\t".
			$ptrGeneData->{$selectedGeneID}{SfirstVariance}."\t".
			$ptrGeneData->{$selectedGeneID}{Dvariance}."\n";
			
		print $filehandle $line;			
	}	
}

sub print_reference_set_of_genes_with_T_for_the_set{
	my ($ptrSelectedGenes, $ptrGeneData, $filehandle)=@_;
	
#	print $filehandle "# S -- number of segregating sites within one spieces in the gene\n";
#	print $filehandle "# D -- number of divergent sites between spieces in the gene\n";
	my $header="";
	my $line="";
	
	$header="#geneID"."\t".
			"feature"."\t".
			"chromosome"."\t".
			"start_position"."\t".
			"end_position"."\t".
			"effective_length_of_gene"."\t".
			"strand"."\t".
			"average_coverage"."\t".
			"theta"."\t".
			"S_first"."\t".
			"S_first_expectation"."\t".
			"S_first_variance"."\t".
			"K"."\t".
			"D"."\t".
			"T"."\n";
		
	print $filehandle $header;
	
	foreach my $selectedGeneID (@{$ptrSelectedGenes}){
		$line= 
			$selectedGeneID."\t". 
			$ptrGeneData->{$selectedGeneID}{feature}."\t".
			$ptrGeneData->{$selectedGeneID}{chromosome}."\t".
			$ptrGeneData->{$selectedGeneID}{start}."\t".
			$ptrGeneData->{$selectedGeneID}{end}."\t".
			$ptrGeneData->{$selectedGeneID}{effectiveLengthOfGene}."\t".
			$ptrGeneData->{$selectedGeneID}{strand}."\t".
			$ptrGeneData->{$selectedGeneID}{avgCoverage}."\t".
			$ptrGeneData->{$selectedGeneID}{theta}."\t".
			$ptrGeneData->{$selectedGeneID}{Sfirst}."\t".
			$ptrGeneData->{$selectedGeneID}{SfirstExpectation}."\t".
			$ptrGeneData->{$selectedGeneID}{SfirstVariance}."\t".
			$ptrGeneData->{$selectedGeneID}{K}."\t".
			$ptrGeneData->{$selectedGeneID}{D}."\t".
			$ptrGeneData->{$selectedGeneID}{T};
		
		$line=$line."\n";
		print $filehandle $line;			
	}	
}

sub load_codon_table{
	my($codon_table_file)=@_;
	my $ptrCodonTable={};

	if (defined $codon_table_file){
		open codonFileHandle,"<",$codon_table_file or die "Could not open codon table file $codon_table_file";
		while (my $line = <codonFileHandle>){
			chomp($line);
			next unless $line;
			next if $line=~ m/^#/;
			my ($codon, $aa) = split /\s*:\s*/, $line;
			$codon=~s/\s+//g;
            $aa=~s/\s+//;
            $codon = uc($codon);
            $aa=uc($aa);
			$ptrCodonTable->{$codon}=$aa;
		}
		close codonFileHandle;
	}else{
		##use default codon table
		$ptrCodonTable = {
			"ATT" => "I", "ATC" => "I", "ATA" => "I", 
			"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L", "TTA" => "L", "TTG" => "L", 
			"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V", 
			"TTT" => "F", "TTC" => "F",
			"ATG" => "M",
			"TGT" => "C", "TGC" => "C",
			"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
			"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G",
			"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P", 
			"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
			"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S", "AGT" => "S", "AGC" => "S",
			"TAT" => "Y", "TAC" => "Y",
			"TGG" => "W",
			"CAA" => "Q", "CAG" => "Q",
			"AAT" => "N", "AAC" => "N",
			"CAT" => "H", "CAC" => "H",
			"GAA" => "E", "GAG" => "E",
			"GAT" => "D", "GAC" => "D",
			"AAA" => "K", "AAG" => "K",
			"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R", "AGA" => "R", "AGG" => "R",
			"TAA" => "-", "TAG" => "-", "TGA" => "-"
		};
	}
	return $ptrCodonTable;
}


1;