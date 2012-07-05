#! /usr/bin/perl

use warnings;
use strict;
use Test::More tests=>113;

use Data::Dumper;

use FindBin qw($RealBin);
use lib "$RealBin/../";

# PopGenTools
use TinaPerl qw(
	_load_line_and_split
	_one_line_string_output
	_calculate_realLength_nextFirstPos_deleteWindowIfNeeded
	_push_data_into_window
	_load_window
	_add_empty_SNP
	_calculate_average_for_window
	average_variance_sliding_window_sync
);

BEGIN{use_ok('VarMath');}
BEGIN{use_ok('VarianceExactCorrection');}
BEGIN{use_ok('TinaPerl');}

#
# 	_load_line_and_split
#

can_ok('TinaPerl', '_load_line_and_split');

my $inString;
my $inFileHandle;
my $chromosome;
my $position;
my $realWindowSize;
my $ptrMeasures=[];
my $ptrAverages=[];
my $outLineString="";
my $ptrWindow=[];
my $nextFirstPosition;

my $firstPosition;
my $firstChromosome; 
my $WINDOW_SIZE = 3;
my $STEP_SIZE = 1;
my $MIN_LENGTH_FRACTION = 0.6;

$inString = "chromosome\tposition\tvar1\tvar2\tvar3\n";
 
open $inFileHandle, "<", \$inString;
($chromosome, $position, $ptrMeasures)=_load_line_and_split($inFileHandle);
close $inFileHandle;

is($chromosome, "chromosome", "chromosome");
is($position, "position", "position");
is(scalar @{$ptrMeasures}, 3, "number of populations");
is($ptrMeasures->[0], "var1", "variance 1");
is($ptrMeasures->[1], "var2", "variance 2");
is($ptrMeasures->[2], "var3", "variance 3");


#
#	_one_line_string_output
#

can_ok('TinaPerl', '_one_line_string_output');

$chromosome = "chromosome";
$position="position";
$realWindowSize = "realWinSize";
$ptrAverages = [ 
				{val => 1, na_ratio => "na_ratio1"},
				{val => 2, na_ratio => "na_ratio2"}
			   ];

$outLineString = _one_line_string_output($chromosome, $position, $realWindowSize, $ptrAverages);

is($outLineString, "chromosome\tposition\trealWinSize\t1\t2\n", "output string");

#
#	_calculate_realLength_nextFirstPos_deleteWindowIfNeeded
#

can_ok('TinaPerl', '_calculate_realLength_nextFirstPos_deleteWindowIfNeeded');


$ptrWindow = [
	{variances => [1,2,3], position=> 1},
	{variances => [1,2,3], position=> 2},
	{variances => [1,2,3], position=> 3},
];

$firstPosition=1;
$position=4;
$firstChromosome = "chromosome";
$chromosome = "chromosome";

($ptrWindow, 
 $realWindowSize, 
 $nextFirstPosition) = _calculate_realLength_nextFirstPos_deleteWindowIfNeeded($ptrWindow,
 																			   $position, 
 																			   $firstPosition, 
 																			   $chromosome,
 																			   $firstChromosome,
 																			   $WINDOW_SIZE,
 																			   $STEP_SIZE,
 																			   $MIN_LENGTH_FRACTION);

is(scalar @{$ptrWindow},3, "window");
is($realWindowSize, 3, "real window size");
is($nextFirstPosition, $firstPosition+$STEP_SIZE, "next starting position of a window");

$ptrWindow = [
	{variances => [1,2,3], position=> 1},
	{variances => [1,2,3], position=> 2},
#	{variances => [1,2,3], position=> 3},
];

$firstPosition=1;
$position=4;
$firstChromosome = "chromosome1";
$chromosome = "chromosome2";

($ptrWindow, 
 $realWindowSize, 
 $nextFirstPosition) = _calculate_realLength_nextFirstPos_deleteWindowIfNeeded($ptrWindow,
 																			   $position, 
 																			   $firstPosition,
 																			   $chromosome,
 																			   $firstChromosome,
 																			   $WINDOW_SIZE,
 																			   $STEP_SIZE,
 																			   $MIN_LENGTH_FRACTION);


is(scalar @{$ptrWindow},2, "window");
is($realWindowSize, 2, "real window size");
is($nextFirstPosition, 2, "next starting position of a window");



$ptrWindow = [
	{variances => [1,2,3], position=> 1},
#	{variances => [1,2,3], position=> 2},
#	{variances => [1,2,3], position=> 3},
];

$firstPosition=1;
$position=4;
$firstChromosome = "chromosome1";
$chromosome = "chromosome2";

($ptrWindow, 
 $realWindowSize, 
 $nextFirstPosition) = _calculate_realLength_nextFirstPos_deleteWindowIfNeeded($ptrWindow,
 																			   $position, 
 																			   $firstPosition, 
 																			   $chromosome,
 																			   $firstChromosome,
 																			   $WINDOW_SIZE,
 																			   $STEP_SIZE,
 																			   $MIN_LENGTH_FRACTION);


is(scalar @{$ptrWindow},0, "window");
is($realWindowSize, 0, "real window size");
is($nextFirstPosition, $position, "next starting position of a window");


#
# _push_data_into_window
#

can_ok('TinaPerl', '_push_data_into_window');

$ptrWindow=[];

$chromosome="chromosome";
$position=1;
$ptrMeasures=["var1","var2", "var3"];
$firstPosition=1;
$WINDOW_SIZE = 4;


$inString = "chromosome\t3\tvar1\tvar2\tvar3\n".
			"chromosome\t4\tvar1\tvar2\tvar3\n".
			"chromosome\t5\tvar1\tvar2\tvar3\n";
open $inFileHandle, "<", \$inString;

($ptrWindow, $chromosome, $position, $ptrMeasures) = _push_data_into_window($ptrWindow, $chromosome, $position, $ptrMeasures, $firstPosition, $inFileHandle, $WINDOW_SIZE);

close $inFileHandle;


is(scalar @{$ptrWindow},3, "window");
is($ptrWindow->[0]{position}, 1, "first SNP");
is($ptrWindow->[1]{position}, 3, "second SNP");
is($ptrWindow->[2]{position}, 4, "third SNP");
is($chromosome, "chromosome", "chromosome");
is($position, 5, "next position");
is(scalar @{$ptrMeasures}, 3, "next measures");
is($ptrMeasures->[0], "var1", "var 1");
is($ptrMeasures->[1], "var2", "var 2");
is($ptrMeasures->[2], "var3", "var 3");

$ptrWindow=[];

$chromosome="chromosome";
$position=2;
$ptrMeasures=["var1","var2", "var3"];
$firstPosition=2;
$WINDOW_SIZE = 3;

$inString = "chromosome\t3\tvar1\tvar2\tvar3\n".
			"chromosome2\t4\tvar1\tvar2\tvar3\n";
open $inFileHandle, "<", \$inString;

($ptrWindow, $chromosome, $position, $ptrMeasures) = _push_data_into_window($ptrWindow, $chromosome, $position, $ptrMeasures, $firstPosition, $inFileHandle, $WINDOW_SIZE);

close $inFileHandle;


is(scalar @{$ptrWindow},2, "window");
is($ptrWindow->[0]{position}, 2, "first SNP");
is($ptrWindow->[1]{position}, 3, "second SNP");
is($chromosome, "chromosome2", "chromosome");
is($position, 4, "next position");
is(scalar @{$ptrMeasures}, 3, "next measures");
is($ptrMeasures->[0], "var1", "var 1");
is($ptrMeasures->[1], "var2", "var 2");
is($ptrMeasures->[2], "var3", "var 3");

#
# _add_empty_SNP
#

can_ok('TinaPerl', '_add_empty_SNP');

$ptrWindow=[];

$ptrWindow=_add_empty_SNP($ptrWindow, 1, 1);

is(scalar @{$ptrWindow},1, "empty array");
is($ptrWindow->[0]{position},1,"first position");
is(scalar @{$ptrWindow->[0]{variances}}, 1, "populations");
is($ptrWindow->[0]{variances}[0], 0, "zero variance");

$ptrWindow=[{variances=>[0.3], position=>3}];

$ptrWindow=_add_empty_SNP($ptrWindow, 1, 1);

is(scalar @{$ptrWindow},2, "non-empty array");
is($ptrWindow->[0]{position},1,"first position");
is(scalar @{$ptrWindow->[0]{variances}}, 1, "populations");
is($ptrWindow->[0]{variances}[0], 0, "zero variance");

#
# _load_window
#

can_ok('TinaPerl', '_load_window');


$ptrWindow=[];

$chromosome="chromosome";
$position=1;
$ptrMeasures=["var1","var2", "var3"];
$firstPosition=1;
$firstChromosome = "chromosome";
$WINDOW_SIZE = 3;
$MIN_LENGTH_FRACTION=0.6;
$STEP_SIZE=2;

$inString = "chromosome\t3\tvar1\tvar2\tvar3\n".
			"chromosome\t4\tvar1\tvar2\tvar3\n".
			"chromosome\t5\tvar1\tvar2\tvar3\n";
			
open $inFileHandle, "<", \$inString;

($ptrWindow, $realWindowSize, $firstPosition, $chromosome, $position, $ptrMeasures)=
_load_window(
$inFileHandle, 
$ptrWindow, 
$chromosome, 
$position, 
$ptrMeasures, 
$firstPosition, 
$firstChromosome, 
$WINDOW_SIZE, 
$MIN_LENGTH_FRACTION, 
$STEP_SIZE);

is($chromosome, "chromosome", "chromosome");
is($firstPosition, 3, "next first position");
is($position, 4, "next position");
is($realWindowSize, 3, "real win size");
is(scalar @{$ptrMeasures}, 3, "measures");
is($ptrMeasures->[0], "var1", "var 1");
is($ptrMeasures->[1], "var2", "var 2");
is($ptrMeasures->[2], "var3", "var 3");
is(scalar @{$ptrWindow},2,"window");
is($ptrWindow->[0]{position}, 1, "first position");
is(scalar @{$ptrWindow->[0]{variances}},3,"first variances");
is($ptrWindow->[0]{variances}[0], "var1", " variance 1");
is($ptrWindow->[0]{variances}[1], "var2", " variance 2");
is($ptrWindow->[0]{variances}[2], "var3", " variance 3");
is($ptrWindow->[1]{position}, 3, "second position");
is($ptrWindow->[1]{variances}[0], "var1", " variance 1");
is($ptrWindow->[1]{variances}[1], "var2", " variance 2");
is($ptrWindow->[1]{variances}[2], "var3", " variance 3");


($ptrWindow, $realWindowSize, $firstPosition, $chromosome, $position, $ptrMeasures)=
_load_window(
$inFileHandle, 
$ptrWindow, 
$chromosome, 
$position, 
$ptrMeasures, 
$firstPosition, 
$firstChromosome, 
$WINDOW_SIZE, 
$MIN_LENGTH_FRACTION, 
$STEP_SIZE);

is($chromosome, undef, "next chromosome");
is($firstPosition, undef, "next first position");
is($realWindowSize, 3, "real win size");
is(scalar @{$ptrMeasures}, 0, "next measures");
is(scalar @{$ptrWindow},3,"window");
is($ptrWindow->[0]{position}, 3, "first position");
is(scalar @{$ptrWindow->[0]{variances}},3,"first variances");
is($ptrWindow->[0]{variances}[0], "var1", " variance 1");
is($ptrWindow->[0]{variances}[1], "var2", " variance 2");
is($ptrWindow->[0]{variances}[2], "var3", " variance 3");
is($ptrWindow->[1]{position}, 4, "second position");
is($ptrWindow->[1]{variances}[0], "var1", " variance 1");
is($ptrWindow->[1]{variances}[1], "var2", " variance 2");
is($ptrWindow->[1]{variances}[2], "var3", " variance 3");
is($ptrWindow->[2]{position}, 5, "third position");
is($ptrWindow->[2]{variances}[0], "var1", " variance 1");
is($ptrWindow->[2]{variances}[1], "var2", " variance 2");
is($ptrWindow->[2]{variances}[2], "var3", " variance 3");

close $inFileHandle;


#
# _calculate_average_for_window
#

my $MIN_COV_FRACTION;
my $ALL_POP_NONZERO;


can_ok('TinaPerl', '_calculate_average_for_window');

$MIN_COV_FRACTION = 1;
$ALL_POP_NONZERO = 1;
$realWindowSize = 3;

$ptrWindow=[
	{variances=>[1,1,1,1], position=>1},
	{variances=>[2,2,2,2], position=>2},
	{variances=>[3,3,3,3], position=>3}
];

$ptrAverages = _calculate_average_for_window($ptrWindow, $realWindowSize, $MIN_COV_FRACTION);

is(scalar @{$ptrAverages}, 4, "number of averages");
is($ptrAverages->[0]{na_ratio}, 0, "na_ratio 1");
is($ptrAverages->[0]{val}, 2, "value 1");
is($ptrAverages->[1]{na_ratio}, 0, "na_ratio 2");
is($ptrAverages->[1]{val}, 2, "value 2");
is($ptrAverages->[2]{na_ratio}, 0, "na_ratio 3");
is($ptrAverages->[2]{val}, 2, "value 3");
is($ptrAverages->[3]{na_ratio}, 0, "na_ratio 4");
is($ptrAverages->[3]{val}, 2, "value 4");


$MIN_COV_FRACTION = 0.6;
$ALL_POP_NONZERO = 0;
$realWindowSize = 3;

$ptrWindow=[
	{variances=>[1,1,1,"na",0], position=>1},
	{variances=>[2,0,"na",2,1], position=>2},
	{variances=>[3,0,3,3,2], position=>3}
];

$ptrAverages = _calculate_average_for_window($ptrWindow, $realWindowSize, $MIN_COV_FRACTION);

is(scalar @{$ptrAverages}, 5, "number of averages");
is($ptrAverages->[0]{na_ratio}, 0, "na_ratio 1");
is($ptrAverages->[0]{val}, 2, "value 1");

is($ptrAverages->[1]{na_ratio}, 0, "na_ratio 2");
is($ptrAverages->[1]{val}, 1/3, "value 2");

is($ptrAverages->[2]{na_ratio}, 1/3, "na_ratio 3");
is($ptrAverages->[2]{val}, 2, "value 3");

is($ptrAverages->[3]{na_ratio}, 1/3, "na_ratio 4");
is($ptrAverages->[3]{val}, 2.5, "value 4");

is($ptrAverages->[4]{na_ratio}, 0, "na_ratio 5");
is($ptrAverages->[4]{val}, 1, "value 5");


$MIN_COV_FRACTION = 0.6;
$ALL_POP_NONZERO = 1;
$realWindowSize = 3;

$ptrWindow=[
	{variances=>[1,1,1,"na",0], position=>1},
	{variances=>[2,0,2,"na",1], position=>2},
	{variances=>[3,0,"na",3,2], position=>3}
];

$ptrAverages = _calculate_average_for_window($ptrWindow, $realWindowSize, $MIN_COV_FRACTION);

is(scalar @{$ptrAverages}, 0, "na_fraction is too high");

$MIN_COV_FRACTION = 0.6;
$ALL_POP_NONZERO = 1;
$realWindowSize = 3;

$ptrWindow=[
	{variances=>[0], position=>1},
	{variances=>[0], position=>2},
	{variances=>[0], position=>3}
];

$ptrAverages = _calculate_average_for_window($ptrWindow, $realWindowSize, $MIN_COV_FRACTION);

is(scalar @{$ptrAverages}, 1, "all populations zero");

is($ptrAverages->[0]{na_ratio}, 0, "na_ratio");
is($ptrAverages->[0]{val}, 0, "value");


#
# average_variance_sliding_window_sync
#


#$inString = "chromosome\t3\t1\t2\t3\n".
#			"chromosome\t4\t1\t2\t3\n".
#			"chromosome\t5\t1\t2\t3\n";

#my $OUTPUT_FILE = "test";

#$WINDOW_SIZE = 3;
#$STEP_SIZE = 2;
#$MIN_COV_FRACTION =1;
#$MIN_LENGTH_FRACTION = 1;
#$ALL_POP_NONZERO = 1;

#average_variance_sliding_window_sync(\$inString, $OUTPUT_FILE, $WINDOW_SIZE, $STEP_SIZE, $MIN_COV_FRACTION, $MIN_LENGTH_FRACTION, $ALL_POP_NONZERO);






