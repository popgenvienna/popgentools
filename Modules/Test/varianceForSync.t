#!/usr/bin/perl

#testing package varianceForSync;

use warnings;
use strict;

use Data::Dumper;
use Test::More tests=>44;
#use Test::More qw(no_plan);

use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/..";

use VarianceExactCorrection;

use TinaPerl qw(
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
);

BEGIN{use_ok('TinaPerl');}
BEGIN{use_ok('VarianceExactCorrection');}


my $MIN_COUNT = 4;
my $MIN_COV = 10;
my $MAX_COV=100;
my $POOL_SIZE=100;

#
#	_population_SNP_counts
#

can_ok('TinaPerl','_population_SNP_counts');

my $ptrH;
my $out;
$ptrH = {A=>4,T=>1,C=>1,G=>1};
$out = _population_SNP_counts($ptrH,$MIN_COUNT);
is($out, 0, "_population_SNP_counts");

$ptrH={A=>4,T=>1,C=>400,G=>1};
$out = _population_SNP_counts($ptrH, $MIN_COUNT);
is($out, 1, "_population_SNP_counts");

$ptrH={A=>0,T=>30,C=>400,G=>4};
$out = _population_SNP_counts($ptrH, $MIN_COUNT);
is($out, 1, "_population_SNP_counts");

$ptrH={A=>4,T=>30,C=>400,G=>4};
$out = _population_SNP_counts($ptrH, $MIN_COUNT);
is($out, 1, "_population_SNP_counts");

#
# 	_population_SNP_coverage
#
can_ok('TinaPerl','_population_SNP_coverage');


$ptrH = {eucov=>20, asterisk=>0};
$out =_population_SNP_coverage($ptrH, $MIN_COV, $MAX_COV);
is($out, 1, "_population_SNP_coverage");

$ptrH = {eucov=>3, asterisk=>0};
$out =_population_SNP_coverage($ptrH, $MIN_COV, $MAX_COV);
is($out, 0, "_population_SNP_coverage");

$ptrH = {eucov=>101, asterisk=>0};
$out =_population_SNP_coverage($ptrH, $MIN_COV, $MAX_COV);
is($out, 0, "_population_SNP_coverage");

$ptrH = {eucov=>22, asterisk=>1};
$out =_population_SNP_coverage($ptrH, $MIN_COV, $MAX_COV);
is($out, 0, "_population_SNP_coverage");

#
# 	_set_na_one_measure_one_population
#
can_ok('TinaPerl', '_set_na_one_measure_one_population');

my $ptrData={};
$ptrH = {chr=>"2L", pos=>1};
my $measure;
$measure="theta";

_set_na_one_measure_one_population($ptrH, $measure, $ptrData);
is(scalar @{$ptrData->{"2L"}[1]{$measure}}, 1,"_set_na_one_measure_one_population");
is($ptrData->{"2L"}[1]{$measure}[0], "na","_set_na_one_measure_one_population");

#
#	_set_na_all_measures_one_population
#
can_ok('TinaPerl', '_set_na_all_measures_one_population');
$ptrData={};
$ptrH = {chr=>"2L", pos=>1};
my $ptrSettings;
$ptrSettings = [{measure=>"pi"}, {measure=>"theta"}, {measure=>"D"}];

_set_na_all_measures_one_population($ptrH, $ptrSettings, $ptrData);
is(scalar @{$ptrData->{"2L"}[1]{"pi"}}, 1,"_set_na_all_measures_one_population");
is(scalar @{$ptrData->{"2L"}[1]{"theta"}}, 1,"_set_na_all_measures_one_population");
is(scalar @{$ptrData->{"2L"}[1]{D}}, 1,"_set_na_all_measures_one_population");
is($ptrData->{"2L"}[1]{"pi"}[0], "na","_set_na_all_measures_one_population");
is($ptrData->{"2L"}[1]{"theta"}[0], "na","_set_na_all_measures_one_population");
is($ptrData->{"2L"}[1]{"D"}[0], "na","_set_na_all_measures_one_population");

#
# 	_set_zero_one_measure_one_population
#
can_ok('TinaPerl','_set_zero_one_measure_one_population');

$ptrData={};
$ptrH = {chr=>"2L", pos=>1};
$measure="theta";

_set_zero_one_measure_one_population($ptrH, $measure, $ptrData);
is(scalar @{$ptrData->{"2L"}[1]{$measure}}, 1,"_set_zero_one_measure_one_population");
is($ptrData->{"2L"}[1]{$measure}[0], 0,"_set_zero_one_measure_one_population");

#
#	_set_zero_all_measures_one_population
#
can_ok('TinaPerl','_set_zero_all_measures_one_population');
$ptrData={};
$ptrH = {chr=>"2L", pos=>1};
$ptrSettings = [{measure=>"pi"}, {measure=>"theta"}, {measure=>"D"}];

_set_zero_all_measures_one_population($ptrH, $ptrSettings, $ptrData);
is(scalar @{$ptrData->{"2L"}[1]{"pi"}}, 1,"_set_zero_all_measures_one_population");
is(scalar @{$ptrData->{"2L"}[1]{"theta"}}, 1,"_set_zero_all_measures_one_population");
is(scalar @{$ptrData->{"2L"}[1]{D}}, 1,"_set_zero_all_measures_one_population");
is($ptrData->{"2L"}[1]{"pi"}[0], 0,"_set_zero_all_measures_one_population");
is($ptrData->{"2L"}[1]{"theta"}[0], 0,"_set_zero_all_measures_one_population");
is($ptrData->{"2L"}[1]{"D"}[0], 0,"_set_zero_all_measures_one_population");


#
#	_calculate_one_measure_one_population
#
can_ok('TinaPerl','_calculate_one_measure_one_population');
$ptrH={A=>10,T=>10,C=>10,G=>10,N=>0,asterisk=>0,eucov=>40, chr=>"2L", pos=>1};
$measure="theta";
$ptrData={};
my $vec = VarianceExactCorrection->new($POOL_SIZE,$MIN_COUNT);
_calculate_one_measure_one_population($ptrH,$measure,$vec,$ptrData);
is($ptrData->{"2L"}[1]{"theta"}[0], 0.427125303468737,"_calculate_one_measure_one_population");

#
#	_calculate_all_measures_one_population
#
can_ok('TinaPerl','_calculate_all_measures_one_population');
$ptrData={};
$ptrSettings = [{measure=>"pi"}, {measure=>"theta"}, {measure=>"D"}];
$ptrH={A=>10,T=>10,C=>10,G=>10,N=>0,asterisk=>0,eucov=>40, chr=>"2L", pos=>1};
$vec = VarianceExactCorrection->new($POOL_SIZE,$MIN_COUNT);

_calculate_all_measures_one_population($ptrH, $ptrSettings, $vec, $ptrData);
is($ptrData->{"2L"}[1]{"theta"}[0], 0.427125303468737,"_calculate_all_measures_one_population");
is($ptrData->{"2L"}[1]{"pi"}[0], 0.90908942402591,"_calculate_all_measures_one_population");
is($ptrData->{"2L"}[1]{"D"}[0], 2.47954271036852,"_calculate_all_measures_one_population");


#
#	_calculate_all_measures_all_populations
#
can_ok('TinaPerl','_calculate_all_measures_all_populations');
$ptrData={};
my $ptrAllPop=[];
$ptrAllPop=[
{A=>10,T=>10,C=>10,G=>10,N=>0,asterisk=>0,eucov=>40, chr=>"2L", pos=>1},
]; 
_calculate_all_measures_all_populations($ptrAllPop, $vec, $ptrSettings, $ptrData, $MIN_COUNT, $MIN_COV, $MAX_COV);
is($ptrData->{"2L"}[1]{"theta"}[0], 0.427125303468737,"_calculate_all_measures_all_populations");
is($ptrData->{"2L"}[1]{"pi"}[0], 0.90908942402591,"_calculate_all_measures_all_populations");
is($ptrData->{"2L"}[1]{"D"}[0], 2.47954271036852,"_calculate_all_measures_all_populations");

#print Dumper($ptrData);
can_ok('TinaPerl','_print_all_data_all_measures');
can_ok('TinaPerl','_print_all_data_one_measure');


#########################
# synchronized_measures #
#########################

my $inFile = "2L\t1\ta\t-\t-\n".
			 "2L\t2\ta\t1:4:4:2:0:0\t1:6:6:6:0:0\n";

open my $filehandle, ">-";
my $ptrOut=[
	{filehandle=>$filehandle,measure=>"pi"},
	{filehandle=>$filehandle,measure=>"theta"},
	{filehandle=>$filehandle,measure=>"D"}
];
synchronized_measures(\$inFile, $ptrOut, $POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV);
close $filehandle;








