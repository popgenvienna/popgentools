#!/usr/bin/perl

package TinaPerl;

use warnings;
use strict;
use Data::Dumper;

#use FindBin qw($RealBin);
#use lib "$RealBin";

use POSIX qw(ceil floor);
use File::Basename;

#export
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
	synchronized_measures
	synchronized_measures_low_memory_requirement
	average_variance_sliding_window_sync
	read_calculate_print_windows_and_log
	read_data_candidates_list_calculate_print_windows_log
);
our @EXPORT_OK = qw(
	_population_SNP_counts
	_population_SNP_coverage
	
	_set_na_one_measure_one_population
	_set_na_all_measures_one_population
	
	_set_zero_one_measure_one_population
	_set_zero_all_measures_one_population
	
	_calculate_one_measure_one_population
	_calculate_all_measures_one_population
	_calculate_all_measures_all_populations
	
	_print_all_data_all_measures
	_print_all_data_one_measure
	
	
	_load_line_and_split
	_one_line_string_output
	_calculate_realLength_nextFirstPos_deleteWindowIfNeeded
	_push_data_into_window
	_load_window
	_add_empty_SNP
	_calculate_average_for_window


	_parse_populations_one_replica_different_conditions
	_load_SNPs
	_load_c
	_create_window_for_selected_SNP
	_process_SNPs
	_print_window
	_print_log_withdrawn_SNPs
);

#PopGenTools
use VarianceExactCorrection;

####################################
#	                               #
# variances for sync pileup file   #
#								   #	
####################################

sub synchronized_measures_low_memory_requirement{
	my ($IN_FILE, $ptrSettings, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV) = @_;
	# $ptrSettings -- array of outputs, where output is a reference to a hash with keys 2: filehandle and measure
	my $vec=VarianceExactCorrection->new($POOL_SIZE,$MIN_COUNT);
	
		open SYNCpileup, "<", $IN_FILE or die "Could not open file $IN_FILE";
	while (my $line = <SYNCpileup>){	
		chomp($line);
		#skip comments
		next if $line =~ m/^#/;

		#decide, weather the position is a SNP;
		my ($chromosome, $position, $reference, @populations) = split "\t", $line;
						
		my $ptrPopulationsSeparately=[];
		
		foreach my $population (@populations){
				my $ptrPopulation={A=>0, T=>0, C=>0, G=>0, N=>0, asterisk=>0, eucov=>0, chr=>$chromosome, pos=>$position};
				
				#skip empty population data
				unless ($population =~ m/\-/){
				
					my @counts = split ":", $population;
				
					$ptrPopulation ={
						A=>$counts[0],
						T=>$counts[1],
						C=>$counts[2],
						G=>$counts[3],
						N=>$counts[4],
						asterisk=>$counts[5],				
						eucov=>$counts[0]+$counts[1]+$counts[2]+$counts[3], 
						chr=>$chromosome,
						pos=>$position
					};
				
				};
				
				push @{$ptrPopulationsSeparately}, $ptrPopulation;
		}
	
		_calculate_and_print_out_all_measures_all_populations($ptrPopulationsSeparately, $vec, $ptrSettings, $MIN_COUNT, $MIN_COV, $MAX_COV);
			#my ($ptrPopulationsSeparately, $vec, $ptrSettings, $MIN_COUNT, $MIN_COV, $MAX_COV)=@_;
	}
	close SYNCpileup;	
}

sub synchronized_measures{
	my ($IN_FILE, $ptrSettings, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV) = @_;
		
	# $ptrSettings -- array of outputs, where output is a reference to a hash with keys 2: filehandle and measure 
		
	my $vec=VarianceExactCorrection->new($POOL_SIZE,$MIN_COUNT);
	
	my $ptrData={};
	
	# $ptrData -- structure for storage all calculated values for different measures 
	# description:
	# $ptrData->{chromosome}[position_on_chromosome]{measure_pi_or_theta_or_D} 
	#            == ptr to array of all measure values, one value per population in an input sync-pileup file 
	
	open SYNCpileup, "<", $IN_FILE or die "Could not open file $IN_FILE";
	while (my $line = <SYNCpileup>){	
		chomp($line);
		#skip comments
		next if $line =~ m/^#/;

		#decide, weather the position is a SNP;
		my ($chromosome, $position, $reference, @populations) = split "\t", $line;
						
		my $ptrPopulationsSeparately=[];
		
		foreach my $population (@populations){
				my $ptrPopulation={A=>0, T=>0, C=>0, G=>0, N=>0, asterisk=>0, eucov=>0, chr=>$chromosome, pos=>$position};
				
				#skip empty population data
				unless ($population =~ m/\-/){
				
					my @counts = split ":", $population;
				
					$ptrPopulation ={
						A=>$counts[0],
						T=>$counts[1],
						C=>$counts[2],
						G=>$counts[3],
						N=>$counts[4],
						asterisk=>$counts[5],				
						eucov=>$counts[0]+$counts[1]+$counts[2]+$counts[3], 
						chr=>$chromosome,
						pos=>$position
					};
				
				};
				
				push @{$ptrPopulationsSeparately}, $ptrPopulation;
		}
	
		_calculate_all_measures_all_populations($ptrPopulationsSeparately, $vec, $ptrSettings, $ptrData, $MIN_COUNT, $MIN_COV, $MAX_COV); 
	}
	close SYNCpileup;
	
	_print_all_data_all_measures($ptrData, $ptrSettings);
}

sub _print_all_data_all_measures{
	my ($ptrData, $ptrSettings)=@_;
	
	foreach my $ptrSetting (@$ptrSettings){
		my $measure = $ptrSetting->{measure};
		my $filehandle = $ptrSetting->{filehandle};
		_print_all_data_one_measure($ptrData, $measure, $filehandle);
	}
}
sub _print_all_data_one_measure{
	my ($ptrData, $measure, $fileHandle)=@_;
	
	print $fileHandle "#calculated measure: $measure\n";
	print $fileHandle "#chromosome\tposition\tmeasures, one for each population\n";
	foreach my $chromosome (keys %$ptrData){
		for (my $position=0; $position< scalar(@{$ptrData->{$chromosome}}); $position++){
			next unless (defined($ptrData->{$chromosome}[$position]));
			
			my $values="";
			my $end = scalar(@{$ptrData->{$chromosome}[$position]{$measure}})-1;
			for (my $i=0; $i<$end; $i++){
				$values=$values.$ptrData->{$chromosome}[$position]{$measure}[$i]."\t";
			}
			$values=$values.$ptrData->{$chromosome}[$position]{$measure}[$end];
			
			print $fileHandle $chromosome,"\t",$position,"\t",$values,"\n"; 
		}
	}
}
#sub _calculate_one_measure_all_populations{
#	my ($ptrPopulationsSeparately,$measure,$vec, $ptrData)=@_;

#	foreach my $ptrPopulation (@{$ptrPopulationsSeparately}){
#		_calculate_one_measure_one_population($ptrPopulation,$measure,$vec,$ptrData);
#	}	
#}

#TODO: this
sub _calculate_and_print_out_all_measures_all_populations{
	my ($ptrPopulationsSeparately, $vec, $ptrSettings, $MIN_COUNT, $MIN_COV, $MAX_COV)=@_;
	
	my $ptrStrings = {};
	
	#initialize values of strings
	foreach my $ptrSetting (@$ptrSettings){
		my $measure = $ptrSetting->{measure};
		$ptrStrings->{$measure}="";
	}
	
	foreach my $ptrPopulation (@$ptrPopulationsSeparately){
		#print Dumper($ptrPopulation);
		if (_population_SNP_coverage($ptrPopulation,$MIN_COV,$MAX_COV)){
			if (_population_SNP_counts($ptrPopulation, $MIN_COUNT)){
				_calculate_and_print_to_strings_all_measures_one_population($ptrPopulation, $ptrSettings, $vec, $ptrStrings);		
				#	my ($ptrPopulation, $ptrSettings, $vec, $ptrStrings)=@_;
			}else{
				_set_zero_and_print_to_strings_all_measures_one_population($ptrPopulation, $ptrSettings, $ptrStrings);	
			}
		}else{
			_set_na_and_print_to_strings_all_measures_one_population($ptrPopulation, $ptrSettings, $ptrStrings);
		}
		#print Dumper($ptrStrings);
	}

	#print all strings to files
	foreach my $ptrSetting (@$ptrSettings){
		my $measure = $ptrSetting->{measure}; 
		my $filehandle = $ptrSetting->{filehandle};
		print $filehandle $ptrStrings->{$measure}."\n";
	}
}

sub _calculate_all_measures_all_populations{
	my ($ptrPopulationsSeparately, $vec, $ptrSettings, $ptrData, $MIN_COUNT, $MIN_COV, $MAX_COV)=@_;
	
	foreach my $ptrPopulation (@$ptrPopulationsSeparately){
		#my $eucov = $ptrPopulation->{eucov};
		if (_population_SNP_coverage($ptrPopulation,$MIN_COV,$MAX_COV)){
			if (_population_SNP_counts($ptrPopulation, $MIN_COUNT)){
				_calculate_all_measures_one_population($ptrPopulation, $ptrSettings, $vec, $ptrData);
			}else{
				_set_zero_all_measures_one_population($ptrPopulation, $ptrSettings, $ptrData);
			}
		}else{
			_set_na_all_measures_one_population($ptrPopulation, $ptrSettings, $ptrData);
		}
	}
}
sub _calculate_all_measures_one_population{
	my ($ptrPopulation, $ptrSettings, $vec, $ptrData)=@_;
	
	foreach my $ptrSetting (@$ptrSettings){
		my $measure = $ptrSetting->{measure};
		_calculate_one_measure_one_population($ptrPopulation, $measure, $vec, $ptrData);
	}	
}


sub _set_na_all_measures_one_population{
	my ($ptrPopulation, $ptrSettings, $ptrData)=@_;
	
	foreach my $ptrSetting (@$ptrSettings){
		#print Dumper($ptrSetting);
		my $measure = $ptrSetting->{measure};
		_set_na_one_measure_one_population($ptrPopulation, $measure, $ptrData);
	}
}
sub _set_na_one_measure_one_population{
	my ($ptrPopulation, $measure, $ptrData)=@_;
	
	my $chromosome = $ptrPopulation->{chr};
	my $position = $ptrPopulation->{pos};
	
	push @{$ptrData->{$chromosome}[$position]{$measure}}, "na";
}

sub _set_na_and_print_to_strings_all_measures_one_population{
	my ($ptrPopulation, $ptrSettings, $ptrStrings)=@_;
	
	foreach my $ptrSetting (@$ptrSettings){
		my $measure = $ptrSetting->{measure};
		my $str = $ptrStrings->{$measure};
		$ptrStrings->{$measure} = _set_na_and_print_to_string_one_measure_one_population($ptrPopulation, $measure, $str);	
	}
	
}

sub _set_na_and_print_to_string_one_measure_one_population{
	my ($ptrPopulation, $measure, $prevStr)=@_;
	
	my $str = $prevStr;
	
	if ($str eq ""){
		my $chromosome = $ptrPopulation->{chr};
		my $position = $ptrPopulation->{pos};
		$str = $chromosome."\t".$position."\tna";
	}else{
		$str=$str."\tna";	
	}
	return $str;
}

sub _set_zero_all_measures_one_population{
	my ($ptrPopulation, $ptrSettings, $ptrData)=@_;

	foreach my $ptrSetting (@$ptrSettings){
		my $measure = $ptrSetting->{measure};
		_set_zero_one_measure_one_population($ptrPopulation, $measure, $ptrData);
	}
}

sub _set_zero_and_print_to_strings_all_measures_one_population{
	my ($ptrPopulation, $ptrSettings, $ptrStrings)=@_;
	
	foreach my $ptrSetting (@$ptrSettings){
		my $measure =  $ptrSetting->{measure};
		my $str = $ptrStrings->{$measure};
		$ptrStrings->{$measure} = _set_zero_and_print_to_string_one_measure_one_population($ptrPopulation, $measure, $str);
	}
	
}

sub _set_zero_and_print_to_string_one_measure_one_population{
	my ($ptrPopulation, $measure, $prevStr)=@_;

	my $str = $prevStr;
	
	if ($str eq ""){
		my $chromosome = $ptrPopulation->{chr};
		my $position = $ptrPopulation->{pos};
		$str = $chromosome."\t".$position."\t0";
	}else{
		$str=$str."\t0";
	}
	return $str;
}

sub _set_zero_one_measure_one_population{
	my ($ptrPopulation, $measure, $ptrData)=@_;
	
	my $chromosome = $ptrPopulation->{chr};
	my $position = $ptrPopulation->{pos};

	push @{$ptrData->{$chromosome}[$position]{$measure}}, 0;	
}
sub _calculate_one_measure_one_population{
	my ($ptrPopulation, $measure, $vec, $ptrData)=@_;
	
	my $value = $vec->calculate_measure($measure, [$ptrPopulation],1);
	
	my $chromosome = $ptrPopulation->{chr};
	my $position = $ptrPopulation->{pos};
	
	push @{$ptrData->{$chromosome}[$position]{$measure}}, $value;
}



sub _calculate_and_print_to_strings_all_measures_one_population{
	my ($ptrPopulation, $ptrSettings, $vec, $ptrStrings)=@_;
	
	foreach my $ptrSetting (@$ptrSettings){
		my $measure =  $ptrSetting->{measure};
		my $str = $ptrStrings->{$measure};
		#print Dumper($ptrStrings, $str);
		$ptrStrings->{$measure} = _calculate_and_print_to_string_one_measure_one_population($ptrPopulation, $measure, $vec, $str);		
	}
	
}

sub _calculate_and_print_to_string_one_measure_one_population{
	my ($ptrPopulation, $measure, $vec, $prevStr)=@_;
	
	my $value = $vec->calculate_measure($measure, [$ptrPopulation],1);
	
	my $str = $prevStr;
	
	#print Dumper($str);
	
	if($str eq ""){
		my $chromosome = $ptrPopulation->{chr};
		my $position = $ptrPopulation->{pos};
		$str= $chromosome."\t".$position."\t".$value;
	}else{
		$str=$str."\t".$value;
	}
	return $str;
}

sub _population_SNP_counts{
	my ($ptrPopulation, $MIN_COUNT) = @_;
	my $tmp=0;
	
	my $a = $ptrPopulation->{A};
	my $t = $ptrPopulation->{T};
	my $c = $ptrPopulation->{C};
	my $g = $ptrPopulation->{G};	
		
	if ($a>=$MIN_COUNT){ $tmp+=1; }
	if ($t>=$MIN_COUNT){ $tmp+=1; }	
	if ($c>=$MIN_COUNT){ $tmp+=1; }
	if ($g>=$MIN_COUNT){ $tmp+=1; }
	
	if ($tmp>=2){
		return 1;
	}else{
		return 0;
	}
}
sub _population_SNP_coverage{
	my ($ptrPopulation, $MIN_COV, $MAX_COV) = @_;
	
	my $eucov = $ptrPopulation->{eucov};
	my $as = $ptrPopulation->{asterisk};
	
	if (($MIN_COV<=$eucov)and($eucov<=$MAX_COV)and($as==0)){
		return 1;	
	}else{
		return 0;
	}
}

#######################################################################
#																	  #				
# average variances for sliding window on synchronized variances file #
#																	  #			
#######################################################################

sub _load_line_and_split{
	my ($InFileHandle) = @_;
	
	my $line;
	my $chromosome;
	my $position;
	my @measures;
	
	if (!eof($InFileHandle)){
		$line = <$InFileHandle>;
    
		chomp($line);
		if ($line ne ""){	
			($chromosome, $position, @measures) = split "\t", $line;
		}else{
			#empty line at the end of file
			$chromosome = undef; $position = undef; @measures = ();
		}
	
		return ($chromosome, $position, \@measures);
	}else{
		$chromosome = undef; $position = undef; @measures = ();
		return ($chromosome, $position, \@measures);
	}
}

sub _one_line_string_output{
	my ($chromosome, $start, $realWindowSize, $ptrAverages, $PRINT_ZERO_LINES, $ALL_POP_NONZERO)=@_;
	
	my $ptrAvgValues = [];
	my $lastPop = scalar @{$ptrAverages} -1;
	for (my $i=0; $i<=$lastPop; $i++){
		push @{$ptrAvgValues}, $ptrAverages->[$i]{val};
	} 
	
	my $ptrNAfractions = [];
	for (my $i=0; $i<=$lastPop; $i++){
		push @{$ptrNAfractions}, $ptrAverages->[$i]{na_ratio};
	}
	
	my $avgValString = join("\t", @{$ptrAvgValues});
	my $naString = join("\t", @{$ptrNAfractions});

	my $zero = 0;
	foreach my $pop (@{$ptrAverages}){			
		if($pop->{val} == 0){
			$zero+=1;
		}
	}
	  	
	my $outString; 	
	
	if (($zero == scalar@{$ptrAverages}) and !$PRINT_ZERO_LINES){
		$outString = "";
	}elsif( ($zero > 0) and ($ALL_POP_NONZERO) ){
		$outString = "";
	}else{
		$outString = $chromosome."\t".$start."\t".$realWindowSize."\t".$avgValString."\n"; #.$naString."\n";
	}
	
	return $outString;
}

sub _calculate_realLength_nextFirstPos_deleteWindowIfNeeded{
	my ($ptrWindow, $position, $firstPosition, $chromosome, $firstChromosome, $WINDOW_SIZE, $STEP_SIZE, $MIN_LENGTH_FRACTION)=@_;
	
	my $actLength;
	my $nextFirstPos;
	
	my $lastWindowIndex = scalar @{$ptrWindow} - 1;
	my $lastPosition = $ptrWindow->[$lastWindowIndex]{position};
	
	# detect why the previous while loop stopped	
	if (! defined($position)){	
		#readed emply line at the end
		if ($lastPosition - $firstPosition + 1< $WINDOW_SIZE*$MIN_LENGTH_FRACTION){
		# too short window, discard it
			$actLength=0;
			$nextFirstPos = $position;
			$ptrWindow=[];
			return ($ptrWindow, $actLength, $nextFirstPos);	
		}else{
		# the window is long enough to calculate averages
			$actLength = $lastPosition - $firstPosition +1;
		
			if ($WINDOW_SIZE-$STEP_SIZE > $WINDOW_SIZE*$MIN_LENGTH_FRACTION){
				$nextFirstPos = $firstPosition + $STEP_SIZE;
			}else{
				$nextFirstPos = $position;
			}		
			return ($ptrWindow, $actLength, $nextFirstPos);
		}
	}elsif(($firstChromosome eq $chromosome)and($position - $firstPosition + 1 > $WINDOW_SIZE)){
		$actLength = $WINDOW_SIZE;
		$nextFirstPos = $firstPosition + $STEP_SIZE;
		return ($ptrWindow, $actLength, $nextFirstPos);
	}elsif(($firstChromosome ne $chromosome)and($lastPosition - $firstPosition + 1< $WINDOW_SIZE*$MIN_LENGTH_FRACTION)){
		$actLength=0;
		$nextFirstPos = $position;
		$ptrWindow=[];
		return ($ptrWindow, $actLength, $nextFirstPos);
	}elsif(($firstChromosome ne $chromosome)and($lastPosition - $firstPosition + 1>= $WINDOW_SIZE*$MIN_LENGTH_FRACTION)){
		$actLength = $lastPosition - $firstPosition +1;
		if ($STEP_SIZE<$WINDOW_SIZE){
			$nextFirstPos = $firstPosition + $STEP_SIZE;
		}else{
			$nextFirstPos = $position;
		}
		return ($ptrWindow, $actLength, $nextFirstPos);
	}	
}

sub _push_data_into_window{
	my ($ptrWindow, $chromosome, $position, $ptrMeasures, $firstPosition, $InFileHandle, $WINDOW_SIZE)=@_;
	
	my $ptrNextMeasures=[];
	my $nextPosition;
	my $nextChromosome;
	my $prevPosition;

	push @{$ptrWindow}, {variances => $ptrMeasures, position => $position};
		
	$prevPosition = $position; 
	($nextChromosome, $nextPosition, $ptrNextMeasures) = _load_line_and_split($InFileHandle);
	
	while ( (defined($nextPosition))and($nextPosition-$firstPosition+1 <= $WINDOW_SIZE)and($chromosome eq $nextChromosome)){
		
		#store data
		push @{$ptrWindow}, {variances => $ptrNextMeasures, position => $nextPosition};	
			
		#read new line
		($nextChromosome, $nextPosition, $ptrNextMeasures) = _load_line_and_split($InFileHandle);
		
		#print Dumper($nextChromosome, $nextPosition, $ptrNextMeasures);
	}
	
	return ($ptrWindow, $nextChromosome, $nextPosition, $ptrNextMeasures);
}


sub _add_empty_SNP{
	my ($ptrWindow, $popNumber, $firstPosition)=@_;
	
	# create empty pop measure array
	my $popLast = $popNumber-1;	
	my $ptrNullMeasures=[];
	for( my $i=0; $i<=$popLast; $i++){
		$ptrNullMeasures->[$i]=0;
	}

	#add first position at the beginning of the array
	unshift(@{$ptrWindow}, ({variances=>$ptrNullMeasures , position=>$firstPosition}));	
	
	return $ptrWindow;
}


sub _load_window{
	my ($InFileHandle, $ptrWindow, $chromosome, $position, $ptrMeasures, $firstPosition, $firstChromosome, $WINDOW_SIZE, $MIN_LENGTH_FRACTION, $STEP_SIZE)=@_;

	my $nextChromosome=$chromosome;
	my $nextPosition = $position;
	my $ptrNextMeasures = $ptrMeasures;
	
	# shift first positions out of the array if array is not empty
	if (scalar @{$ptrWindow} != 0){
		while ((scalar @{$ptrWindow} > 0)and($ptrWindow->[0]{position} < $firstPosition)){
			shift(@{$ptrWindow});
		}
	}
		
	#add first position if it is not alredy there
	if ( ((scalar @{$ptrWindow} == 0)and($firstPosition != $position))
	   or((scalar @{$ptrWindow} != 0)and($firstPosition != $ptrWindow->[0]{position})) 
	){
		$ptrWindow=_add_empty_SNP($ptrWindow, scalar @{$ptrMeasures}, $firstPosition);	
	}
	
	
	my $realWindowSize;
	my $nextFirstPos;
	
	#if position from lastly read line is in the WINDOW_SIZE distance on the same chromosome add the position and variances to the array
	if (($position-$firstPosition+1 <= $WINDOW_SIZE)and($chromosome eq $firstChromosome)){
		
		($ptrWindow, $nextChromosome, $nextPosition, $ptrNextMeasures)= _push_data_into_window($ptrWindow, $chromosome, $position, $ptrMeasures, $firstPosition, $InFileHandle, $WINDOW_SIZE);
				
		#print Dumper($ptrWindow, $nextPosition, $firstPosition, $nextChromosome, $chromosome);
		
		
		($ptrWindow, 
		 $realWindowSize,
		 $nextFirstPos) = _calculate_realLength_nextFirstPos_deleteWindowIfNeeded( $ptrWindow, 
		 																		   $nextPosition, 
		 																		   $firstPosition, 
		 																	 	   $nextChromosome, 
		 																	       $chromosome,
		 																	 	   $WINDOW_SIZE,
		 																	       $STEP_SIZE,
		 																	       $MIN_LENGTH_FRACTION);
		 																	       
		return ($ptrWindow, $realWindowSize, $nextFirstPos, $nextChromosome, $nextPosition, $ptrNextMeasures);	
		
	}elsif(($position-$firstPosition+1 > $WINDOW_SIZE)and($chromosome eq $firstChromosome)){
	# the same chromosome, jumt to position first + stepsize 	
		$realWindowSize = $ptrWindow->[scalar @{$ptrWindow} - 1]{position} - $firstPosition +1;
		$nextFirstPos = $firstPosition + $STEP_SIZE;
		return ($ptrWindow, $realWindowSize, $nextFirstPos, $nextChromosome, $nextPosition, $ptrNextMeasures);
	}else{
	#chromosome switch	
		$ptrWindow = [];
		$realWindowSize = 0;
		$nextFirstPos = 1;	
		return ($ptrWindow, $realWindowSize, $nextFirstPos, $nextChromosome, $nextPosition, $ptrNextMeasures)		
	}
}


sub _calculate_average_for_window{
	my ($ptrWindow, $realWindowSize, $MIN_COV_FRACTION)=@_;
	
	my $ptrAverages =[];
	my $ptrSum=[];
	
	if ($realWindowSize == 0){
		return $ptrAverages;
		
	}else{

		my $last = scalar @{$ptrWindow} -1;
		my $lastPop = scalar@{$ptrWindow->[0]{variances}}-1; 

		#initialize sums to 0
		for(my $j=0; $j<=$lastPop; $j++){
			$ptrSum->[$j]={val => 0, na_count => 0};
		}
	
	
		#calculate sum of variances in window
		for (my $j=0; $j<=$lastPop; $j++){
			for (my $i=0; $i<=$last; $i++){
				if ( $ptrWindow->[$i]{variances}[$j] ne "na"){
					$ptrSum->[$j]{val}+=$ptrWindow->[$i]{variances}[$j];
				}else{
					$ptrSum->[$j]{na_count}+=1;
				}
			}
		}
		
		#check min_cov_fraction for each population
		my $continue = 1;
	
		foreach my $ptrPopulation (@{$ptrSum}){
			if ($ptrPopulation->{na_count} / $realWindowSize > (1-$MIN_COV_FRACTION) ){
				$continue = 0;
			}
		}
	
		#setting average values
		if ($continue){
			for (my $j=0; $j<=$lastPop; $j++){
				if ($realWindowSize == $ptrSum->[$j]{na_count}){
					$ptrAverages->[$j]{val} = 0;
				}else{
					$ptrAverages->[$j]{val} = $ptrSum->[$j]{val} / ($realWindowSize - $ptrSum->[$j]{na_count});
				}
				$ptrAverages->[$j]{na_ratio} = $ptrSum->[$j]{na_count} / $realWindowSize;
			}
		}
	
		return $ptrAverages;
	}
}



sub average_variance_sliding_window_sync{
	my ($SYNC_VARIANCE_FILE, $OUTPUT_FILE, $WINDOW_SIZE, $STEP_SIZE, $MIN_COV_FRACTION, $MIN_LENGTH_FRACTION, $PRINT_ZERO_LINES, $ALL_POP_NONZERO)=@_;
	
	open my $inFileHandle, "<", $SYNC_VARIANCE_FILE or die "Could not open sync-variance file $SYNC_VARIANCE_FILE";
	open my $outFileHandle, ">", $OUTPUT_FILE or die "Could not open output file $OUTPUT_FILE";


	my $firstPosition=1;
	my $chromosome;
	my $position; 
	my $ptrMeasures = [];
	my $ptrWindow = [];

	if (!eof($inFileHandle)){ 
		($chromosome, $position, $ptrMeasures)=_load_line_and_split($inFileHandle);

		while (!eof($inFileHandle)){

			#print $position, "\n";
			
			my $realWindowSize;
			my $actChromosome = $chromosome;
	
			( $ptrWindow, 
	  		  $realWindowSize, 
	  	 	  $firstPosition, 
		  	  $chromosome,
		  	  $position, 
	      	  $ptrMeasures) =  _load_window( $inFileHandle, 
	  										$ptrWindow,
	  										$chromosome,
	  										$position, 
	  										$ptrMeasures, 
	  										$firstPosition, 
	  										$chromosome, 
	  										$WINDOW_SIZE,
		  									$MIN_LENGTH_FRACTION,
		  									$STEP_SIZE);
      		
    		my $ptrAverages = _calculate_average_for_window($ptrWindow, $realWindowSize, $MIN_COV_FRACTION);
	    	
    		if (scalar @{$ptrAverages} != 0){
    			my $outString = _one_line_string_output($actChromosome, $ptrWindow->[0]{position}, $realWindowSize, $ptrAverages, $PRINT_ZERO_LINES, $ALL_POP_NONZERO);
    			print $outFileHandle $outString;
    		}
		}
	}
	close $outFileHandle;
	close $inFileHandle;
}



######################################################
#
# linkage disequilibrium comparison ?? (heloise)
#
######################################################


sub _parse_populations_one_replica_different_conditions{
	my ($earlyCountsString, $lateCountsString, $ptrWithdrawnSNPs, $chromosome, $position, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $IGNORE)=@_;
	
	my @earlyCounts=(0,0,0,0,0,0); 
	my @lateCounts=(0,0,0,0,0,0);
	
	unless ($earlyCountsString =~ m/\-/){
		@earlyCounts = split ":", $earlyCountsString;	
	}
	unless ($lateCountsString =~ m/\-/){
		@lateCounts = split ":", $lateCountsString;	
	}
	 
	
	my $domAllele = "";
	my $recAllele = "";
	my $countDomAllele = 0;
	my $countRecAllele = 0;
	my $countDomEarly = 0;
	my $countRecEarly = 0;
	
	my $condition;
	
	if (defined($IGNORE) ){
		$condition = (($earlyCounts[4]<=$IGNORE)and($earlyCounts[5]<=$IGNORE)and($lateCounts[4]<=$IGNORE)and($lateCounts[5]<=$IGNORE));
	}else{
		$condition = (($earlyCounts[4]==0)and($earlyCounts[5]==0)and($lateCounts[4]==0)and($lateCounts[5]==0));
	}
	
	
	if ($condition){
			
		my %SumCounts=();		
		$SumCounts{"A"}=$earlyCounts[0]+$lateCounts[0];
		$SumCounts{"T"}=$earlyCounts[1]+$lateCounts[1];
		$SumCounts{"C"}=$earlyCounts[2]+$lateCounts[2];
		$SumCounts{"G"}=$earlyCounts[3]+$lateCounts[3];

		my %counts_early=();
		$counts_early{"A"}=$earlyCounts[0];
		$counts_early{"T"}=$earlyCounts[1];
		$counts_early{"C"}=$earlyCounts[2];
		$counts_early{"G"}=$earlyCounts[3];
		
		my $ptrSortedByValue = [];
		foreach my $nucleotide (sort {$SumCounts{$b}<=> $SumCounts{$a}} (keys(%SumCounts))){
			my %hash = (val=> $SumCounts{$nucleotide}, nucleotide => $nucleotide);
			push @{$ptrSortedByValue}, \%hash;
		}
		
		if ($ptrSortedByValue->[2]{val} <= $MAX_COUNT_FOR_3RD_ALLELE){		
			if (
				($ptrSortedByValue->[0]{val}>=$MIN_COUNT_FOR_SELECTED_ALLELES)and
				($ptrSortedByValue->[1]{val}>=$MIN_COUNT_FOR_SELECTED_ALLELES)
			){
				my $nucleotide1 =  $ptrSortedByValue->[0]{nucleotide};
				my $nucleotide2 =  $ptrSortedByValue->[1]{nucleotide};
				if (
					($counts_early{$nucleotide1}>0)and 
					($counts_early{$nucleotide2}>0)and
					($SumCounts{$nucleotide1}-$counts_early{$nucleotide1}>0)and
					($SumCounts{$nucleotide2}-$counts_early{$nucleotide2}>0)
				){
				
				$countDomAllele = $ptrSortedByValue->[0]{val};
				$countRecAllele = $ptrSortedByValue->[1]{val};		
				
				$domAllele= $ptrSortedByValue->[0]{nucleotide};
				$recAllele= $ptrSortedByValue->[1]{nucleotide};
				
				$countDomEarly = $counts_early{$domAllele};
				$countRecEarly = $counts_early{$recAllele};
				}else{
					push @{$ptrWithdrawnSNPs}, 
						$chromosome."\t".
						$position."\t".
						"At least one of selected alleles counts equals to zero in one of the populations.";	
				}
			}else{
				push @{$ptrWithdrawnSNPs}, 
					$chromosome."\t".
					$position."\t".
					"1st allele count: ".$ptrSortedByValue->[0]{val}."\t".
					"2nd allele count: ".$ptrSortedByValue->[1]{val}."\t".
					"MIN_COUNT_FOR_SELECTED_ALLELES: ".$MIN_COUNT_FOR_SELECTED_ALLELES;
			}						
		}else{
			push @{$ptrWithdrawnSNPs}, 
				$chromosome."\t".
				$position."\t".
				"3rd allele count: ".$ptrSortedByValue->[2]{val}."\t".
				"MAX_COUNT_FOR_3RD_ALLELE: ".$MAX_COUNT_FOR_3RD_ALLELE;
		}
	}else{
		push @{$ptrWithdrawnSNPs}, 
			$chromosome."\t".
			$position."\t".
			"MORE than ".$IGNORE." InDels or unkNown nucleotides occured in pileup file on this position.";	
	}	
	return ($domAllele, $recAllele, $countDomAllele, $countRecAllele, $countDomEarly, $countRecEarly, $ptrWithdrawnSNPs);
}

sub _load_SNPs{
	my ($SNP_CANDIDATES_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $IGNORE, $LANE)=@_;
	
	my $ptrAllSNPs = {};
	my $ptrWithdrawnSNPs=[];
	my $ptrAllAbsValues = {};
		
	open snpCandidatesHandle, "<", $SNP_CANDIDATES_FILE;
	while (my $line = <snpCandidatesHandle>){
		chomp($line);
		my (
			$chromosome,  
			$position, 
			$reference, 
			$earlyCounts1, 
			$earlyCounts3, 
			$lateCounts1, 
			$lateCounts3
		) = split "\t", $line;
	
		my $earlyCounts;
		my $lateCounts;
			
		if ($LANE == 1){
			$earlyCounts = $earlyCounts1;
			$lateCounts = $lateCounts1;
		}elsif($LANE == 3){
			$earlyCounts = $earlyCounts3;
			$lateCounts = $lateCounts3;
		}

		my (
			$A_Allele, 
			$a_Allele, 
			$counts_A_together, 
			$counts_a_together,
			$counts_A_early,
			$counts_a_early,
			$ptrWithdrawnSNPs
		) = _parse_populations_one_replica_different_conditions(
				$earlyCounts, 
				$lateCounts, 
				$ptrWithdrawnSNPs, 
				$chromosome, 
				$position,
				$MIN_COUNT_FOR_SELECTED_ALLELES,
				$MAX_COUNT_FOR_3RD_ALLELE, 
				$IGNORE,
			);	
			
			
	
		# if correct SNP found, create a record in hash $ptrAllSNPs 
		if ($counts_A_together > 0){
		
			my $counts_A_late = $counts_A_together - $counts_A_early;
			my $counts_a_late = $counts_a_together - $counts_a_early;
	
			my $diff_in_freq_A = $counts_A_late/($counts_A_late + $counts_a_late) - $counts_A_early/($counts_A_early + $counts_a_early);
			my $diff_in_freq_a = $counts_a_late/($counts_A_late + $counts_a_late) - $counts_a_early/($counts_A_early + $counts_a_early);
	
			my $diff_in_freq_abs_max=0;
			
			$diff_in_freq_abs_max = abs($diff_in_freq_A);
	
			my $overall_freq_A = $counts_A_together/ ($counts_A_together + $counts_a_together);	
			my $overall_freq_a = $counts_a_together/ ($counts_A_together + $counts_a_together);
		
			my $key = $chromosome."_".$position;
			
			my %tmpHash = (
				nucleotide_allele_A =>$A_Allele, 
				nucleotide_allele_a =>$a_Allele,
				diff_in_freq_A => $diff_in_freq_A,
				diff_in_freq_a => $diff_in_freq_a,
#				diff_in_freq_abs_max => $diff_in_freq_abs_max,
				overall_freq_A => $overall_freq_A,
				overall_freq_a => $overall_freq_a,
				counts_late_A => $counts_A_late,
				counts_late_a => $counts_a_late,
				counts_early_A => $counts_A_early,
				counts_early_a => $counts_a_early
			);
			$ptrAllSNPs->{$chromosome}[$position]=\%tmpHash;
			
			$ptrAllAbsValues->{$key}=$diff_in_freq_abs_max;
		}
	}
	close snpCandidatesHandle;
	
	return ($ptrAllSNPs, $ptrAllAbsValues, $ptrWithdrawnSNPs); 
}


sub _load_c{
	my ($C_FILE)=@_;
	my $ptrC = {};
	# window length == 10^5 
	open inFileHandle, "<", $C_FILE;
	
	while (my $line = <inFileHandle>){
		chomp($line);
		next if $line =~ m/^#/;
		my ($chromosome, $start, $c) = split "\t", $line;
		$ptrC->{$chromosome}{$start}=$c*10**(-6);
	}
	
	close inFileHandle;

	return $ptrC;
}

sub _load_candidates{
	my ($CANDIDATES_FILE)=@_;
	
	my $ptrCandidates={};
	
	open inFileHandle, "<", $CANDIDATES_FILE;
	my $i=0;
	while (my $line =<inFileHandle>){
		chomp($line);
		next if $line =~ m/^#/;
		my ($chromosome, $position) = split "\t", $line;
		my $key = $chromosome."_".$position;
		$ptrCandidates->{$key}=$i;
		$i+=1;
	}
	close inFileHandle;
	return $ptrCandidates;
}


sub _create_window_for_selected_SNP{
	my ($ptrAllSNPs, $selectedChromosome, $selectedPosition, $WINDOW_SIZE, $ptrC, $ptrDiscarded)=@_;
	
	my $ptrWindow=[];
	my $ptrDiscardedUp = $ptrDiscarded;
	my $ptrChromosomeSNPs = [];
	my $WINDOW_SIZE_HALF=$WINDOW_SIZE/2;
	
	if (!defined($ptrAllSNPs->{$selectedChromosome}[$selectedPosition])){
		push @$ptrDiscardedUp, "$selectedChromosome\t$selectedPosition\t This position is not defined in sync pileup or had been withdrawn according to script settings.";
	}else{
		my $lastChromosomeSNPposition = (scalar @{$ptrAllSNPs->{$selectedChromosome}}) -1;
		$ptrChromosomeSNPs = $ptrAllSNPs->{$selectedChromosome};
	 
		#print Dumper($WINDOW_SIZE_HALF, $lastChromosomeSNPposition, $ptrChromosomeSNPs);
 	
 		my $first;
 		my $last;
 	
		if ((0<$selectedPosition-$WINDOW_SIZE_HALF+1)and($selectedPosition + $WINDOW_SIZE_HALF <= $lastChromosomeSNPposition)){
			$first = $selectedPosition - $WINDOW_SIZE_HALF+1;
			$last = $selectedPosition + $WINDOW_SIZE_HALF; 
		}elsif((0<$selectedPosition-$WINDOW_SIZE_HALF+1)and($selectedPosition + $WINDOW_SIZE_HALF > $lastChromosomeSNPposition)){
			$first = $selectedPosition - $WINDOW_SIZE_HALF+1;
			$last = $lastChromosomeSNPposition;
		}elsif((0>=$selectedPosition-$WINDOW_SIZE_HALF+1)and($selectedPosition + $WINDOW_SIZE_HALF <= $lastChromosomeSNPposition)){
			$first = 1;
			$last = $selectedPosition + $WINDOW_SIZE_HALF;
		}else{
			$first = 1;
			$last = $lastChromosomeSNPposition;
		}
	
		my $overallFreqSelected = $ptrAllSNPs->{$selectedChromosome}[$selectedPosition]{overall_freq_A};
		my $diffInFreqSelected = $ptrAllSNPs->{$selectedChromosome}[$selectedPosition]{diff_in_freq_A};
		
		my $domSelected = $ptrAllSNPs->{$selectedChromosome}[$selectedPosition]{nucleotide_allele_A};
		
		my $ptrSelected = $ptrAllSNPs->{$selectedChromosome}[$selectedPosition];
		
		for (my $index = $first; $index<= $last; $index++){
			# g=0.01
			# t=400
			# r^2 = 0.02 + 0.337/(1+c*dist+2gt(1-exp(-dist/t)))
			next unless defined($ptrChromosomeSNPs->[$index]);
			my $posInWin = $index - $selectedPosition;
			my $dist = abs($index-$selectedPosition);
			my $key = floor($index/10**5) * 10**5 ;
			my $c = $ptrC->{$selectedChromosome}{$key};	if(!defined($c)){$c=0};
			my $r_sqr = 0.02 + 0.337/(1+$c*$dist + 8*(1-exp(-$dist/400)));
			
			
			my $ptrOther = $ptrAllSNPs->{$selectedChromosome}[$index];
			
			my $domOther = $ptrOther->{nucleotide_allele_A};
			
			my $overallFreqOther;
			my $diffInFreqOther;
			
			#print Dumper($c, $dist, $r_sqr, $domSelected, $domOther);
			
			if ($domSelected eq $domOther){
				 $overallFreqOther = $ptrAllSNPs->{$selectedChromosome}[$index]{overall_freq_A};
				 $diffInFreqOther= $ptrAllSNPs->{$selectedChromosome}[$index]{diff_in_freq_A};
			}else{
				$overallFreqOther = $ptrAllSNPs->{$selectedChromosome}[$index]{overall_freq_a};
				$diffInFreqOther = $ptrAllSNPs->{$selectedChromosome}[$index]{diff_in_freq_a};
			}
		
			#print Dumper($overallFreqSelected, $overallFreqOther, $diffInFreqSelected, $diffInFreqOther);
			
			my $q_AB_min_q_aB = sqrt( $r_sqr * ($overallFreqOther * (1-$overallFreqOther)) / ($overallFreqSelected* (1-$overallFreqSelected)) );
			
			my $expectedDiffInFreq = $diffInFreqSelected * $q_AB_min_q_aB;
			
			my $p_a;
			my $p_b; 
			my $direction;
			if ($ptrSelected->{overall_freq_a} >= $ptrOther->{overall_freq_a}){
				$p_a = $ptrSelected->{overall_freq_a};
				$p_b = $ptrOther->{overall_freq_a};
				$direction = 1; 
			}else{
				$p_a = $ptrOther->{overall_freq_a};
				$p_b = $ptrSelected->{overall_freq_a};
				$direction = -1;
			}
			
			my $upper_bound_for_r_sqr = ($p_b*(1-$p_a)) / ($p_a*(1-$p_b)); 
			
			my $q_AB_min_q_aB_upper_bound = sqrt($upper_bound_for_r_sqr * ($overallFreqOther * (1-$overallFreqOther)) / ($overallFreqSelected* (1-$overallFreqSelected)) );
			
			my $diffInFreqExpUpperBound = $q_AB_min_q_aB_upper_bound * $diffInFreqSelected; 
			
			
			my %tmpHash=(
				chromosome => $selectedChromosome,
				positionOnChromosome => $index,
				positionInWindow => $posInWin,
				r_sqr => $r_sqr,
				upper_bound_for_r_sqr => $upper_bound_for_r_sqr,
				p_a => $p_a,
				p_b => $p_b,
				direction => $direction,
				q_AB_min_q_aB => $q_AB_min_q_aB,
				q_AB_min_q_aB_upper_bound => $q_AB_min_q_aB_upper_bound,
				diff_in_freq_expectation_upper_bound =>$diffInFreqExpUpperBound,
				diff_in_freq => $diffInFreqOther,
				diff_in_freq_expectation => $expectedDiffInFreq,
				overall_freq_selected=>$overallFreqSelected,
				overall_freq => $overallFreqOther,
				diff_in_freq_selected => $diffInFreqSelected,	 
			);
			push @{$ptrWindow}, \%tmpHash;	
		}
	}			
	return ($ptrWindow, $ptrDiscardedUp);
}

sub _print_window{
	my ($ptrWindow, $fileName)=@_;
	
	open outFileHandle, ">", $fileName;
	print outFileHandle "#chromosome\tpositionOnChromosome\tpositionInWindowSelected0\tdiffInFreqSelected\tdiffInFreqReal\tdiffInFreqExpectation\tdiffInFreqExpectationUpperBound\tdiffInDiffs\toverallFreq\toverallFreqSelected\tq_AB_min_q_aB\tq_AB_min_q_aB_upper_bound\tr_sqr\tupper_bound_r_sqr\tp_a\tp_b\tdirection\n";
	
	my $ptrRecord={};
	
	foreach $ptrRecord (@$ptrWindow){
		my $chromosome = $ptrRecord->{chromosome};
		my $position = $ptrRecord->{positionOnChromosome};
		my $posInWin = $ptrRecord->{positionInWindow};
		
		my $diffSelected = $ptrRecord->{diff_in_freq_selected};
		my $diff = $ptrRecord->{diff_in_freq};
		my $diffExp = $ptrRecord->{diff_in_freq_expectation};
		my $diffExpUpperBound = $ptrRecord->{diff_in_freq_expectation_upper_bound};

		my $overallFreq = $ptrRecord->{overall_freq};
		my $overallFreqSelected = $ptrRecord->{overall_freq_selected};
		
		my $q_AB_min_q_aB = $ptrRecord->{q_AB_min_q_aB};
		my $q_AB_min_q_aB_upper_bound = $ptrRecord->{q_AB_min_q_aB_upper_bound};
		
		my $r_sqr = $ptrRecord->{r_sqr};
		my $upper_bound_for_r_sqr = $ptrRecord->{upper_bound_for_r_sqr};

		my $p_a = $ptrRecord->{p_a};
		my $p_b = $ptrRecord->{p_b};
		my $direction = $ptrRecord->{direction};

		my $diffInDiffs;
		
		if (abs($diff) - abs($diffExpUpperBound) > 0){
			$diffInDiffs = abs($diff) - abs($diffExpUpperBound);
		}else{
			$diffInDiffs = 0;
		}

		print outFileHandle $chromosome."\t".
							$position."\t".
							$posInWin."\t".
							$diffSelected."\t".
							$diff."\t".
							$diffExp."\t".
							$diffExpUpperBound."\t".
							$diffInDiffs."\t".
							$overallFreq."\t".
							$overallFreqSelected."\t".
							$q_AB_min_q_aB."\t".
							$q_AB_min_q_aB_upper_bound."\t".
							$r_sqr."\t".
							$upper_bound_for_r_sqr."\t".
							$p_a."\t".
							$p_b."\t".
							$direction."\n";
			
	};
};

sub _print_log_withdrawn_SNPs{
	my ($ptrWithdrawnSNPs, $logFileName, $ptrDiscarded)=@_;
	
	open logFileHandle, ">", $logFileName;
	print logFileHandle "#withdrawn SNPs while loading data from sync pileup:\n";
	print logFileHandle "#chromosome\tposition\treason\n";
	my $record;
	foreach $record (@$ptrWithdrawnSNPs){
		print logFileHandle $record."\n";
	}
	print logFileHandle "#\n";

	print logFileHandle "#not processed selected SNPs:\n";
	print logFileHandle "#chromosome\tposition\treason\n";
	foreach $record (@$ptrDiscarded){
		print logFileHandle $record."\n";
	}
	
	close logFileHandle; 
	
}

sub _process_SNPs{
	my ($ptrAllSNPs, $ptrAllAbsValues, $WINDOW_SIZE, $ptrC, $OUT)=@_;

	my $i=0;
	my $ptrWindow =[];
	my $ptrDiscarded = [];
	my $list = $OUT.".list";
	open listHandle, ">", $list;
	foreach my $selectedSNPidentifier ( sort {$ptrAllAbsValues->{$b}<=>$ptrAllAbsValues->{$a}} keys(%$ptrAllAbsValues) ){
		my ($chromosome, $position) = split "_", $selectedSNPidentifier;
		if($i < 10){			
			($ptrWindow, $ptrDiscarded) = _create_window_for_selected_SNP($ptrAllSNPs, $chromosome, $position, $WINDOW_SIZE, $ptrC, $ptrDiscarded);
			
			my $prevOut = $OUT;
			$OUT = $OUT.$i.".selectedSNP"."-chr-".$chromosome."-pos-".$position."-windowSize-".$WINDOW_SIZE.".window";
			_print_window($ptrWindow, $OUT);
			$OUT = $prevOut;
			$i+=1;
		}	
		print listHandle $i."\t".$chromosome."\t".$position."\n";
	}
	close listHandle;
	return $ptrDiscarded;
}

sub _process_SNPs_from_list{
	my ($ptrAllSNPs, $WINDOW_SIZE, $ptrC, $ptrCandidates, $OUT)=@_;
	
	my $i=0;
	my $ptrWindow =[];
	my $ptrDiscarded = [];
	foreach my $key ( sort {$ptrCandidates->{$a}<=>$ptrCandidates->{$b}} keys(%$ptrCandidates) ){
		my ($chromosome, $position) = split "_", $key;
		($ptrWindow, $ptrDiscarded) = _create_window_for_selected_SNP($ptrAllSNPs, $chromosome, $position, $WINDOW_SIZE, $ptrC, $ptrDiscarded);
		
		my $prevOut = $OUT;
		$OUT = $OUT.$i.".candidateSNP"."-chr-".$chromosome."-pos-".$position."-windowSize-".$WINDOW_SIZE.".window";
		_print_window($ptrWindow, $OUT);
		$OUT = $prevOut;
		$i+=1;
	}
	return $ptrDiscarded;
}

sub read_calculate_print_windows_and_log{
	my ($SNP_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $C_FILE, $OUT, $WINDOW_SIZE, $IGNORE, $LANE)=@_;
	
	my ($ptrAllSNPs, $ptrAllAbsValues, $ptrWithdrawnSNPs) = _load_SNPs($SNP_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $IGNORE, $LANE);
	my $ptrC = _load_c($C_FILE);
	
	_process_SNPs($ptrAllSNPs, $ptrAllAbsValues, $WINDOW_SIZE, $ptrC, $OUT);
	
	$OUT=$OUT.".log";
	_print_log_withdrawn_SNPs($ptrWithdrawnSNPs, $OUT);
}


#sub read_calculate_print_windows_and_log_no_recomb{
#	my ($SNP_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $OUT, $WINDOW_SIZE)=@_;
#	my ($ptrAllSNPs, $ptrAllAbsValues, $ptrWithdrawnSNPs) = _load_SNPs($SNP_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE);
#	
#}

sub read_data_candidates_list_calculate_print_windows_log{
	my ($SNP_FILE, $CANDIDATES_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $C_FILE, $OUT, $WINDOW_SIZE, $IGNORE, $LANE)=@_;

	my ($ptrAllSNPs, $ptrAllAbsValues, $ptrWithdrawnSNPs) = _load_SNPs($SNP_FILE, $MIN_COUNT_FOR_SELECTED_ALLELES, $MAX_COUNT_FOR_3RD_ALLELE, $IGNORE, $LANE);
	my $ptrC = _load_c($C_FILE);
	my $ptrCandidates = _load_candidates($CANDIDATES_FILE);
	
	
	my $ptrDiscarded = [];
	$ptrDiscarded = _process_SNPs_from_list($ptrAllSNPs, $WINDOW_SIZE, $ptrC, $ptrCandidates, $OUT);
	$OUT=$OUT.".log";
	_print_log_withdrawn_SNPs($ptrWithdrawnSNPs, $OUT, $ptrDiscarded);	
}




1;
