#!/usr/bin/perl

use strict;
use warnings; 
use Test::More tests => 299;
#use Test::More qw(no_plan);
use Data::Dumper;
use FindBin qw($RealBin);
use lib "$RealBin/..";

use GffGtfParser qw(
	_correct_chromosome_name
	 
	_check_chromosome_arms 
	_update_chromosome_arms_list 
	
	_load
	_overwrite_feature_whole_annotation 
	
	_add_feature_to_specified_area
	_add_feature_whole_annotation
	_is_in_featuresArray
	_add_features
	
	_overwrite_and_add_feat_hash 
	
	_invert_feature_hash 
	_compute_feature_lengths_GFF 
	_calculate_characteristics_GFF
	
	_add_gene_ID_to_annotation_only_CDS_GTF
	_add_gene_ID_to_hash_only_CDS_GTF
	_add_gene_IDs_only_CDS_GTF
);

BEGIN{use_ok('GffGtfParser');}	


#
# _correct_chromosome_name
#

can_ok('GffGtfParser', '_correct_chromosome_name');


my $chr; 

$chr="chr2L";
$chr =  _correct_chromosome_name($chr);
is($chr, "2L", "_correct_chromosome_name, corrected chromosome name");

$chr="2L";
$chr =  _correct_chromosome_name($chr);
is($chr, "2L", "_correct_chromosome_name, original chromosome name");

#
# _check_chromosome_arms
#

can_ok('GffGtfParser', '_check_chromosome_arms');

my $ignore_intergenic;
my $ptrChromosomeArms={};
my $all_chromosome_arms_ok;
my $ptrChromosomeArmsMissing =[];

$ignore_intergenic = 1;
$all_chromosome_arms_ok = _check_chromosome_arms($ignore_intergenic, $ptrChromosomeArms, $ptrChromosomeArmsMissing);
is($all_chromosome_arms_ok,1,"_check_chromosome_arms, ignore_intergenic");
is(scalar@{$ptrChromosomeArmsMissing},0,"_check_chromosome_arms, ptrChromosomeArmsMissing");

$ignore_intergenic=undef;
$ptrChromosomeArmsMissing =[];
$ptrChromosomeArms = {"2L"=>1, "2R"=>0};
$all_chromosome_arms_ok = _check_chromosome_arms($ignore_intergenic, $ptrChromosomeArms, $ptrChromosomeArmsMissing);
is($all_chromosome_arms_ok,0,"_check_chromosome_arms, do not ignore intergenic");
is(scalar@{$ptrChromosomeArmsMissing},1,"_check_chromosome_arms, ptrChromosomeArmsMissing size");
is($ptrChromosomeArmsMissing->[0],"2R","_check_chromosome_arms, ptrChromosomeArmsMissing content");

$ignore_intergenic=0;
$ptrChromosomeArmsMissing =[];
$ptrChromosomeArms = {"2L"=>0, "2R"=>1};
$all_chromosome_arms_ok = _check_chromosome_arms($ignore_intergenic, $ptrChromosomeArms, $ptrChromosomeArmsMissing);
is($all_chromosome_arms_ok,0,"_check_chromosome_arms, do not ignore intergenic");
is(scalar@{$ptrChromosomeArmsMissing},1,"_check_chromosome_arms, ptrChromosomeArmsMissing size");
is($ptrChromosomeArmsMissing->[0],"2L","_check_chromosome_arms, ptrChromosomeArmsMissing content");

#
# _update_chromosome_arms_list
#

can_ok('GffGtfParser', '_update_chromosome_arms_list');

my $feat;
$ptrChromosomeArms={};
$ignore_intergenic=0; 

$chr="2L";
$feat = "CDS";
_update_chromosome_arms_list($ptrChromosomeArms, $ignore_intergenic, $chr, $feat);
is(scalar(keys %$ptrChromosomeArms), 1, "_update_chromosome_arms_list, ptrChromosomeArms size");
is($ptrChromosomeArms->{"2L"}, 0, "_update_chromosome_arms_list, ptrChromosomeArms value");

$chr="2L";
$feat = "chromosome_arm";
_update_chromosome_arms_list($ptrChromosomeArms, $ignore_intergenic, $chr, $feat);
is(scalar keys %$ptrChromosomeArms, 1, "_update_chromosome_arms_list, ptrChromosomeArms size");
is($ptrChromosomeArms->{"2L"}, 1, "_update_chromosome_arms_list, ptrChromosomeArms value");

$chr="2R";
$feat = "chromosome_arm";
_update_chromosome_arms_list($ptrChromosomeArms, $ignore_intergenic, $chr, $feat);
is(scalar keys %$ptrChromosomeArms, 2, "_update_chromosome_arms_list, ptrChromosomeArms size");
is($ptrChromosomeArms->{"2L"}, 1, "_update_chromosome_arms_list, ptrChromosomeArms value");
is($ptrChromosomeArms->{"2R"}, 1, "_update_chromosome_arms_list, ptrChromosomeArms value");

$chr="2R";
$feat = "intron";
_update_chromosome_arms_list($ptrChromosomeArms, $ignore_intergenic, $chr, $feat);
is(scalar keys %$ptrChromosomeArms, 2, "_update_chromosome_arms_list, ptrChromosomeArms size");
is($ptrChromosomeArms->{"2L"}, 1, "_update_chromosome_arms_list, ptrChromosomeArms value");
is($ptrChromosomeArms->{"2R"}, 1, "_update_chromosome_arms_list, ptrChromosomeArms value");


#
# _load 
#

can_ok('GffGtfParser','_load');

my $gff;
my $rGFFlist;

$gff = "2L\tFlyBase\tchromosome_arm\t1\t3\t.\t.\t.";
$rGFFlist = _load(\$gff, \%GffGtfParser::featHash);
is(scalar(@$rGFFlist),1,"gff line loaded");
is($rGFFlist->[0]{chromosome},"2L", "gff line");
is($rGFFlist->[0]{feat},"chromosome_arm", "gff line");
is($rGFFlist->[0]{start},"1", "gff line");
is($rGFFlist->[0]{end},"3", "gff line");
is($rGFFlist->[0]{score},".", "gff line");
is($rGFFlist->[0]{strand},".", "gff line");
is($rGFFlist->[0]{offset},".", "gff line");


$gff = "#comment line\n".
 	   "2L\tFlyBase\tchromosome_arm\t1\t3\t.\t.\t.";
$rGFFlist = _load(\$gff, \%GffGtfParser::featHash);
is(scalar(@{$rGFFlist}),1,"gff file with a comment line loaded");

my $ptrFeatures;
$ptrFeatures = {"chromosome_arm"=>1};
$gff = "2L\tFlyBase\tchromosome_arm\t1\t3\t.\t.\t.\n".
	   "2L\tFlyBase\tintron\t1\t3\t.\t.\t.";	
$rGFFlist = _load(\$gff, $ptrFeatures);
is(scalar(@{$rGFFlist}),1,"gff file with an irrelevant feature");
is($rGFFlist->[0]{feat},"chromosome_arm", "--- correct feature loaded");


$gff = "2L\tFlyBase\tchromosome_arm\t1\t3\t.\t.\t.\n".
	   "2L\tFlyBase\tintron\t1\t3\t.\t.\t.";
$ptrFeatures = {"chromosome_arm"=>1, "intron"=> 1};
$rGFFlist = _load(\$gff, $ptrFeatures, 0);
is(scalar(@{$rGFFlist}),2,"gff file, ignore_intergenic=0");

$gff = "2L\tFlyBase\tintron\t1\t3\t.\t.\t.";
$ptrFeatures = {"intron"=> "i"};
$rGFFlist = _load(\$gff, $ptrFeatures, 1);
is(scalar(@{$rGFFlist}),1,"gff file, ignore_intergenic=1");
is($rGFFlist->[0]{feat},"intron", "gff file");

my $is_gtf;
my $ptrGeneIDs={};

$is_gtf=1;
my $gtf = "chr4\tdm3_flyBaseGene\tCDS\t257895\t258185\t0.000000\t\+\t0\tgene_id \"CG1674-RB\"; transcript_id \"CG1674-RB\"";
$ptrFeatures = {"CDS"=> "C"};
$rGFFlist = _load(\$gtf, $ptrFeatures, 1, $is_gtf);
is(scalar @$rGFFlist, 1, "gtf line loaded");
is($rGFFlist->[0]{chromosome},"4", "gtf chromosome name");
is($rGFFlist->[0]{feat},"CDS", "gtf feature");
is($rGFFlist->[0]{start},"257895", "gtf feature start");
is($rGFFlist->[0]{end},"258185", "gtf feature end");
is($rGFFlist->[0]{score},"0.000000", "gtf score");
is($rGFFlist->[0]{strand},"+", "gtf strand");
is($rGFFlist->[0]{offset},"0", "gtf offset");
is($rGFFlist->[0]{geneID},"CG1674-RB", "gtf geneID");
is($rGFFlist->[0]{transcriptID},"CG1674-RB", "gtf transcriptID");


#
# _overwrite_feature_whole_annotation 
#

can_ok('GffGtfParser', '_overwrite_feature_whole_annotation');
my $rAnnotation;

my $gff4 = "2L\tFlyBase\tchromosome_arm\t1\t3\t.\t.\t.\n";
my $gff5 = "2L\tFlyBase\tintron\t2\t3\t.\t.\t.\n";
my $gff6 = "2L\tFlyBase\texon\t2\t2\t.\t.\t.\n";

my $gff7 = "2R\tFlyBase\tchromosome_arm\t1\t5\t.\t.\t.\n";
my $gff8 = "2R\tFlyBase\tintron\t1\t3\t.\t.\t.\n";
my $gff9 = "2R\tFlyBase\texon\t2\t2\t.\t.\t.\n";

# 1 chromosome, string $gff4;
$rAnnotation={};
$rGFFlist = _load(\$gff4, \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm", \%GffGtfParser::featHash);
is(scalar(keys %$rAnnotation),1, "number of chromosomes");
my $sum=0;
foreach (@{$rAnnotation->{"2L"}}){
	next unless defined($_);
	if ($_->{feat} eq "*"){
		$sum+=1;
	}
}
is($sum,3,"number of intergenic nucleotides");

# 1 chromosome, string $gff4.$gff5;
$rAnnotation={};
$rGFFlist = _load(\($gff4.$gff5), \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "intron");

is($rAnnotation->{"2L"}[1]{feat}, "*", "intergenic position");
is($rAnnotation->{"2L"}[2]{feat}, "i", "intron position");
is($rAnnotation->{"2L"}[3]{feat}, "i", "intron position");

# 1 chromosome, string $gff4.$gff5.$gff6;
$rAnnotation={};
$rGFFlist = _load(\($gff4.$gff5.$gff6), \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "intron");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "exon");

is($rAnnotation->{"2L"}[1]{feat}, "*", "intergenic position");
is($rAnnotation->{"2L"}[2]{feat}, "e", "exon position");
is($rAnnotation->{"2L"}[3]{feat}, "i", "intron position");

# 2 chromosomes, string $gff4.$gff7;
$rAnnotation={};
$rGFFlist = _load(\($gff4.$gff7), \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm");

is(scalar(keys %$rAnnotation), 2, "number of chromosomes, 2 chromosomes: 2L, 2R");
$sum=0;
foreach (@{$rAnnotation->{"2L"}}){
	next unless defined($_);
	if ($_->{feat} eq "*"){
		$sum+=1;
	}
}
is($sum,3,"number of intergenic nucleotides in 2L");
$sum=0;
foreach (@{$rAnnotation->{"2R"}}){
	next unless defined($_);
	if ($_->{feat} eq "*"){
		$sum+=1;
	}
}
is($sum,5,"number of intergenic nucleotides in 2R");

# 2 chromosomes, string $gff4.$gff5.$gff6.$gff7.$gff8.$gff9
$rAnnotation={};
$rGFFlist = _load(\($gff4.$gff5.$gff6.$gff7.$gff8.$gff9), \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "intron");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "exon");

is($rAnnotation->{"2L"}[1]{feat}, "*", "intergenic position, 2L");
is($rAnnotation->{"2L"}[2]{feat}, "e", "exon position, 2L");
is($rAnnotation->{"2L"}[3]{feat}, "i", "intron position, 2L");

is($rAnnotation->{"2R"}[1]{feat}, "i", "intron position, 2R");
is($rAnnotation->{"2R"}[2]{feat}, "e", "exon position, 2R");
is($rAnnotation->{"2R"}[3]{feat}, "i", "intron position, 2R");
is($rAnnotation->{"2R"}[4]{feat}, "*", "intergenic position, 2R");
is($rAnnotation->{"2R"}[5]{feat}, "*", "intergenic position, 2R");



#
# _add_feature_to_specified_area
#

can_ok('GffGtfParser', '_add_feature_to_specified_area');

my $ptrChromosome = [undef, undef, {feat=>"c"}, {feat=>"C"}, {feat=>"*"} ];

my $ptrAnnotation = {"2L"=>$ptrChromosome};

_add_feature_to_specified_area($ptrAnnotation, "2L", 1, 4, "C");

is($ptrAnnotation->{"2L"}[2]{feat}, "cC", "_add_feature_to_specified_area, different code");
is($ptrAnnotation->{"2L"}[3]{feat}, "C", "_add_feature_to_specified_area, the same code");
is($ptrAnnotation->{"2L"}[4]{feat}, "*C", "_add_feature_to_specified_area, non-character code");

#
# _add_feature_whole_annotation
#

can_ok('GffGtfParser', '_add_feature_whole_annotation');

my $gff10 = "2L\tFlyBase\tCDS\t1\t2\t.\t.\t.\n";
my $gff11 = "2L\tFlyBase\tpseudogene\t1\t3\t.\t.\t.\n";
my $gff12 = "2L\tFlyBase\tCDS\t2\t3\t.\t.\t.\n";

my $gff13 = "2R\tFlyBase\tCDS\t1\t2\t.\t.\t.\n";
my $gff14 = "2R\tFlyBase\tpseudogene\t1\t3\t.\t.\t.\n";
my $gff15 = "2R\tFlyBase\tmiRNA\t2\t3\t.\t.\t.\n";
my $gff16 = "2R\tFlyBase\tmiRNA\t4\t5\t.\t.\t.\n";
my $gff17 = "2R\tFlyBase\tmiRNA\t2\t5\t.\t.\t.\n";

# 1 chromosome, string $gff4.$gff5.$gff10;
$rAnnotation={};
$rGFFlist = _load(\($gff4.$gff5.$gff6.$gff10), \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "intron");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "exon");

_add_feature_whole_annotation($rAnnotation, $rGFFlist, "CDS");

is($rAnnotation->{"2L"}[1]{feat}, "*C", "intergenic + CDS, 2L");
is($rAnnotation->{"2L"}[2]{feat}, "eC", "exon + CDS, 2L");
is($rAnnotation->{"2L"}[3]{feat}, "i", "intron, 2L");

# 1 chromosome, string $gff4.$gff5.$gff6.$gff10.$gff11.$gff12;
$rAnnotation={};
$rGFFlist = _load(\($gff4.$gff5.$gff6.$gff10.$gff11.$gff12), \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "intron");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "exon");

_add_feature_whole_annotation($rAnnotation, $rGFFlist, "CDS");
_add_feature_whole_annotation($rAnnotation, $rGFFlist, "pseudogene");

like($rAnnotation->{"2L"}[1]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/p/, "pseudogene, 2L");

like($rAnnotation->{"2L"}[2]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/p/, "pseudogene, 2L");

like($rAnnotation->{"2L"}[3]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/p/, "pseudogene, 2L");

# 2 chromosomes, string $gff4.$gff5.$gff6.$gff7.$gff8.$gff9.$gff10.$gff11.$gff12.$gff13.$gff14.$gff15.$gff16.$gff17;
$rAnnotation={};
$rGFFlist = _load(\($gff4.$gff5.$gff6.$gff7.$gff8.$gff9.$gff10.$gff11.$gff12.$gff13.$gff14.$gff15.$gff16.$gff17), \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "intron");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "exon");

_add_feature_whole_annotation($rAnnotation, $rGFFlist, "CDS");
_add_feature_whole_annotation($rAnnotation, $rGFFlist, "pseudogene");
_add_feature_whole_annotation($rAnnotation, $rGFFlist, "miRNA");

like($rAnnotation->{"2R"}[1]{feat}, qr/C/, "CDS, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/p/, "pseudogene, 2R");

like($rAnnotation->{"2R"}[2]{feat}, qr/C/, "CDS, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/m/, "miRNA, 2R");

like($rAnnotation->{"2R"}[3]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/p/, "pseudogene, 2R");

like($rAnnotation->{"2R"}[4]{feat}, qr/m/, "miRNA, 2R");

like($rAnnotation->{"2R"}[5]{feat}, qr/m/, "miRNA, 2R");

#
# _is_in_featuresArray
#

can_ok('GffGtfParser', '_is_in_featuresArray');

my $feature="CDS"; 
my $ptrFeaturesArray;

$ptrFeaturesArray = ["chromosome_arm", "CDS"];
is(_is_in_featuresArray($feature, $ptrFeaturesArray), 1, "is in features array");

$ptrFeaturesArray = [];
is(_is_in_featuresArray($feature, $ptrFeaturesArray), 0, "is not in features array");

$ptrFeaturesArray = ["intron", "exon"];
is(_is_in_featuresArray($feature, $ptrFeaturesArray), 0, "is not in features array");

#
# _add_features
#

$rAnnotation={};
$rGFFlist = _load(\($gff4.$gff5.$gff6.$gff7.$gff8.$gff9.$gff10.$gff11.$gff12.$gff13.$gff14.$gff15.$gff16.$gff17), \%GffGtfParser::featHash);
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "chromosome_arm");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "intron");
_overwrite_feature_whole_annotation($rAnnotation, $rGFFlist, "exon");

_add_features($rAnnotation, $rGFFlist, ["CDS",  "pseudogene", "miRNA"]);

like($rAnnotation->{"2L"}[1]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/p/, "pseudogene, 2L");

like($rAnnotation->{"2L"}[2]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/p/, "pseudogene, 2L");

like($rAnnotation->{"2L"}[3]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/p/, "pseudogene, 2L");

like($rAnnotation->{"2R"}[1]{feat}, qr/C/, "CDS, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/p/, "pseudogene, 2R");

like($rAnnotation->{"2R"}[2]{feat}, qr/C/, "CDS, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/m/, "miRNA, 2R");

like($rAnnotation->{"2R"}[3]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/p/, "pseudogene, 2R");

like($rAnnotation->{"2R"}[4]{feat}, qr/m/, "miRNA, 2R");

like($rAnnotation->{"2R"}[5]{feat}, qr/m/, "miRNA, 2R");


my $gffAllFeatures = 
                  "2L\tFlyBase\tchromosome_arm\t1\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tintron\t2\t3\t.\t.\t.\n".
				  "2L\tFlyBase\texon\t2\t2\t.\t.\t.\n".
				  "2L\tFlyBase\tncRNA\t1\t1\t.\t.\t.\n".
				  "2L\tFlyBase\tCDS\t1\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tCDS\t2\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tfive_prime_UTR\t1\t2\t.\t.\t.\n".
				  "2L\tFlyBase\tthree_prime_UTR\t2\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tenhancer\t3\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tfive_prime_UTR\t1\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tpre_miRNA\t1\t1\t.\t.\t.\n".				  
				  "2L\tFlyBase\tpseudogene\t3\t3\t.\t.\t.\n".				  				  				  
				  "2R\tFlyBase\tchromosome_arm\t1\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tintron\t1\t3\t.\t.\t.\n".
				  "2R\tFlyBase\ttRNA\t2\t2\t.\t.\t.\n".
				  "2R\tFlyBase\tsnoRNA\t5\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tsnRNA\t1\t1\t.\t.\t.\n".
				  "2R\tFlyBase\trRNA\t4\t4\t.\t.\t.\n".
				  "2R\tFlyBase\tthree_prime_UTR\t3\t3\t.\t.\t.\n".
				  "2R\tFlyBase\tfive_prime_UTR\t2\t4\t.\t.\t.\n".				  
				  "2R\tFlyBase\tmiRNA\t1\t3\t.\t.\t.\n".
				  "2R\tFlyBase\tenhancer\t2\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tregulatory_region\t3\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tregulatory_region\t4\t4\t.\t.\t.\n".
				  "2R\tFlyBase\tpseudogene\t1\t3\t.\t.\t.\n".
				  "2R\tFlyBase\ttransposable_element\t2\t3\t.\t.\t.\n".
				  "2R\tFlyBase\tpseudogene\t2\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tpre_miRNA\t1\t1\t.\t.\t.\n".
				  "2R\tFlyBase\tpre_miRNA\t4\t5\t.\t.\t.\n".
  				  "2R\tFlyBase\tmiRNA\t5\t5\t.\t.\t.\n";

# 2L: 
# 1: 1CFa
# 2: eCFT
# 3: iCThFp
#
# 2R:
# 1: 4mpa
# 2: 2Fmhpt
# 3: iTFmhrpt
# 4: 5Fhrpa
# 5: 3mhrpa



#
# _overwrite_and_add_feat_hash
#

my $ptrGff;

can_ok('GffGtfParser', '_overwrite_and_add_feat_hash');

$rAnnotation={};
$ptrGff=[];
	
$ptrGff = _load(\$gffAllFeatures, \%GffGtfParser::featHash);

#print Dumper($ptrGff);

_overwrite_and_add_feat_hash($rAnnotation, $ptrGff);

is(scalar(keys %$rAnnotation),2, "number of chromosomes");

like($rAnnotation->{"2L"}[1]{feat}, qr/1/, "ncRNA, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/a/, "pre_miRNA, 2L");

like($rAnnotation->{"2L"}[2]{feat}, qr/e/, "exon, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/T/, "three_prime_UTR, 2L");

like($rAnnotation->{"2L"}[3]{feat}, qr/i/, "intron, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/T/, "three_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/h/, "three_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/p/, "pseudogene, 2L");

like($rAnnotation->{"2R"}[1]{feat}, qr/4/, "snRNA, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/a/, "pre_miRNA, 2R");

like($rAnnotation->{"2R"}[2]{feat}, qr/2/, "tRNA, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/t/, "transposable_elenemt, 2R");

like($rAnnotation->{"2R"}[3]{feat}, qr/i/, "intron, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/T/, "three_prime_UTR, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/t/, "transposable_elenemt, 2R");

like($rAnnotation->{"2R"}[4]{feat}, qr/5/, "rRNA, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/a/, "pre_miRNA, 2R");

like($rAnnotation->{"2R"}[5]{feat}, qr/3/, "snoRNA, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/a/, "pre_miRNA, 2R");

#
# _invert_feature_hash 
#

can_ok('GffGtfParser', '_invert_feature_hash');

my $rInverse = {};
$rInverse = _invert_feature_hash();

my $inverseOK = 1;
my $inverseKeys=0;
foreach my $code (keys %$rInverse){
	my $tmpFeat = $rInverse->{$code};
	$inverseKeys+=1;
	
	if ($code eq "*"){
		next if $rInverse->{$code} eq "intergenic";
		$inverseOK=0;
	}else{
		next if $code eq $GffGtfParser::featHash{$tmpFeat};
		$inverseOK=0;
	}
}
is($inverseOK,1,"inverse hash");
is($inverseKeys, 17, "inverse count");


#
# _compute_feature_lengths 
#

can_ok('GffGtfParser', '_compute_feature_lengths_GFF');

$rAnnotation={};
$rAnnotation= load_GFF_file_with_chromosome_arms(\$gffAllFeatures);

my $rFeatureLengths = _compute_feature_lengths_GFF($rAnnotation, $rInverse);

is($rFeatureLengths->{"intergenic"},0,"intergenic counts");
is($rFeatureLengths->{"intron"},2,"intron counts");
is($rFeatureLengths->{"exon"},1,"counts");
is($rFeatureLengths->{"ncRNA"},1,"counts");
is($rFeatureLengths->{"tRNA"},1,"tRNA counts");
is($rFeatureLengths->{"snoRNA"},1,"snoRNA counts");
is($rFeatureLengths->{"snRNA"},1,"snRNA counts");
is($rFeatureLengths->{"rRNA"},1,"rRNA counts");
is($rFeatureLengths->{"CDS"},3,"CDS counts");
is($rFeatureLengths->{"five_prime_UTR"},6,"five_prime_UTR counts");
is($rFeatureLengths->{"three_prime_UTR"},3,"three_prime_UTR counts");
is($rFeatureLengths->{"enhancer"},5,"enhancer counts");
is($rFeatureLengths->{"miRNA"},4,"mRNA counts");
is($rFeatureLengths->{"regulatory_region"},3,"regulatory_region counts");
is($rFeatureLengths->{"pseudogene"},6,"pseudogene counts");
is($rFeatureLengths->{"transposable_element"},2,"transposable_element counts");
is($rFeatureLengths->{"pre_miRNA"},4,"pre_miRNA counts");

#
# _calculate_characteristics 
#

can_ok('GffGtfParser', '_calculate_characteristics_GFF');

my $PILEUP_FILE = "2L\t1\tN\t8\tCcGgcCcG\taaaaaaaa\n".
				  "2L\t2\tN\t8\tCcGgcCcG\taaaaaaaa\n".
				  "2L\t3\tN\t8\tCcGgcCcG\taaaaaaaa\n".
				  "2R\t1\tN\t8\tCcGgcCcG\taaaaaaaa\n".
				  "2R\t2\tN\t8\tCcGgcCcG\taaaaaaaa\n".
				  "2R\t3\tN\t8\tCcGgcCcG\taaaaaaaa\n".
				  "2R\t4\tN\t8\tCcGgcCcG\taaaaaaaa\n".
				  "2R\t5\tN\t8\tCcGgcCcG\taaaaaaaa\n";
				  
###default settings:
my $QUAL_ENCODING = "illumina";
my $MIN_COUNT = 2;
my $MIN_COV = 2;
my $MAX_COV = 100000;
my $MIN_QUAL = 1;
my $POOL_SIZE=100;

my $rNumbers=_calculate_characteristics_GFF($rAnnotation, 
		\$PILEUP_FILE, 
		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL, 
		$POOL_SIZE);

is($rNumbers->{CDS}{pi}, 0.750105011201128,"pi");
is($rNumbers->{CDS}{theta}, 0.689766164273035,"theta");
is($rNumbers->{CDS}{D}, 0.684940610557322,"D"); 

#
# 	_add_gene_ID_to_annotation_only_CDS_GTF
#
my $geneID;
my $chromosome;
my $start;
my $end;

can_ok('GffGtfParser', '_add_gene_ID_to_annotation_only_CDS_GTF');

#my ($ptrAnnotation, $geneID, $chromosome, $start, $end)=@_;
$ptrAnnotation={};
$geneID="myGeneID";
$chromosome = "2L";
$start=2;
$end=2;

_add_gene_ID_to_annotation_only_CDS_GTF($ptrAnnotation, $geneID, $chromosome, $start, $end);
is($ptrAnnotation->{"2L"}[1], undef, "_add_gene_ID_to_annotation_only_CDS_GTF, undef position");
is($ptrAnnotation->{"2L"}[2]{geneID}[0], $geneID, "_add_gene_ID_to_annotation_only_CDS_GTF, geneID ok");

#
#	_add_gene_ID_to_hash_only_CDS_GTF
#
can_ok('GffGtfParser', '_add_gene_ID_to_hash_only_CDS_GTF');

$ptrGeneIDs ={};
$geneID="myGeneID";
$chromosome = "2L";
$start=2;
$end=2;
_add_gene_ID_to_hash_only_CDS_GTF($ptrGeneIDs, $geneID, $chromosome, $start, $end);
is($ptrGeneIDs->{$geneID}[0]{chromosome}, "2L", "_add_gene_ID_to_hash_only_CDS_GTF, chromosome ok");
is($ptrGeneIDs->{$geneID}[0]{start}, 2, "_add_gene_ID_to_hash_only_CDS_GTF, start pos ok");
is($ptrGeneIDs->{$geneID}[0]{end}, 2, "_add_gene_ID_to_hash_only_CDS_GTF, end pos ok");

#
#	_add_gene_IDs_only_CDS_GTF
#

can_ok('GffGtfParser', '_add_gene_IDs_only_CDS_GTF');


$ptrAnnotation={};
$ptrGeneIDs ={};
$geneID="myGeneID";
$chromosome = "2L";
$start=2;
$end=2;

my $ptrGtf=[{chromosome=>$chromosome, feat=>"CDS", start=>$start, end=>$end, geneID=>$geneID}];

_add_gene_IDs_only_CDS_GTF($ptrAnnotation, $ptrGtf, $ptrGeneIDs);
#print Dumper($ptrAnnotation, $ptrGeneIDs);
is($ptrAnnotation->{"2L"}[1], undef, "_add_gene_ID_only_CDS_GTF, undef position");
is($ptrAnnotation->{"2L"}[2]{geneID}[0], $geneID, "_add_gene_ID_only_CDS_GTF, geneID ok");

is($ptrGeneIDs->{$geneID}[0]{chromosome}, "2L", "_add_gene_ID_only_CDS_GTF, chromosome ok");
is($ptrGeneIDs->{$geneID}[0]{start}, 2, "_add_gene_ID_only_CDS_GTF, start pos ok");
is($ptrGeneIDs->{$geneID}[0]{end}, 2, "_add_gene_ID_only_CDS_GTF, end pos ok");

#
# load_GFF_file_with_chromosome_arms #
#

can_ok('GffGtfParser', 'load_GFF_file_with_chromosome_arms');

# 2 chromosomes string $gffAllFeatures
$rAnnotation={};
$rAnnotation = load_GFF_file_with_chromosome_arms(\$gffAllFeatures);

is(scalar(keys %$rAnnotation),2, "number of chromosomes");

like($rAnnotation->{"2L"}[1]{feat}, qr/1/, "ncRNA, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/a/, "pre_miRNA, 2L");

like($rAnnotation->{"2L"}[2]{feat}, qr/e/, "exon, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/T/, "three_prime_UTR, 2L");

like($rAnnotation->{"2L"}[3]{feat}, qr/i/, "intron, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/T/, "three_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/h/, "three_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/p/, "pseudogene, 2L");

like($rAnnotation->{"2R"}[1]{feat}, qr/4/, "snRNA, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/a/, "pre_miRNA, 2R");

like($rAnnotation->{"2R"}[2]{feat}, qr/2/, "tRNA, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/t/, "transposable_elenemt, 2R");

like($rAnnotation->{"2R"}[3]{feat}, qr/i/, "intron, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/T/, "three_prime_UTR, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/t/, "transposable_elenemt, 2R");

like($rAnnotation->{"2R"}[4]{feat}, qr/5/, "rRNA, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/a/, "pre_miRNA, 2R");

like($rAnnotation->{"2R"}[5]{feat}, qr/3/, "snoRNA, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/a/, "pre_miRNA, 2R");


my $gffNoArms = 
             #     "2L\tFlyBase\tchromosome_arm\t1\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tintron\t2\t3\t.\t.\t.\n".
				  "2L\tFlyBase\texon\t2\t2\t.\t.\t.\n".
				  "2L\tFlyBase\tncRNA\t1\t1\t.\t.\t.\n".
				  "2L\tFlyBase\tCDS\t1\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tCDS\t2\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tfive_prime_UTR\t1\t2\t.\t.\t.\n".
				  "2L\tFlyBase\tthree_prime_UTR\t2\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tenhancer\t3\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tfive_prime_UTR\t1\t3\t.\t.\t.\n".
				  "2L\tFlyBase\tpre_miRNA\t1\t1\t.\t.\t.\n".				  
				  "2L\tFlyBase\tpseudogene\t3\t3\t.\t.\t.\n".				  				  				  
				#  "2R\tFlyBase\tchromosome_arm\t1\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tintron\t1\t3\t.\t.\t.\n".
				  "2R\tFlyBase\ttRNA\t2\t2\t.\t.\t.\n".
				  "2R\tFlyBase\tsnoRNA\t5\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tsnRNA\t1\t1\t.\t.\t.\n".
				  "2R\tFlyBase\trRNA\t4\t4\t.\t.\t.\n".
				  "2R\tFlyBase\tthree_prime_UTR\t3\t3\t.\t.\t.\n".
				  "2R\tFlyBase\tfive_prime_UTR\t2\t4\t.\t.\t.\n".				  
				  "2R\tFlyBase\tmiRNA\t1\t3\t.\t.\t.\n".
				  "2R\tFlyBase\tenhancer\t2\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tregulatory_region\t3\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tregulatory_region\t4\t4\t.\t.\t.\n".
				  "2R\tFlyBase\tpseudogene\t1\t3\t.\t.\t.\n".
				  "2R\tFlyBase\ttransposable_element\t2\t3\t.\t.\t.\n".
				  "2R\tFlyBase\tpseudogene\t2\t5\t.\t.\t.\n".
				  "2R\tFlyBase\tpre_miRNA\t1\t1\t.\t.\t.\n".
				  "2R\tFlyBase\tpre_miRNA\t4\t5\t.\t.\t.\n".
  				  "2R\tFlyBase\tmiRNA\t5\t5\t.\t.\t.\n";
  				  
#
# load_GFF_file_without_chromosome_arms 
#

can_ok('GffGtfParser', 'load_GFF_file_without_chromosome_arms');

# 2 chromosomes string $gffAllFeatures
$rAnnotation={};
$rAnnotation = load_GFF_file_without_chromosome_arms(\$gffNoArms);

is(scalar(keys %$rAnnotation),2, "number of chromosomes");

like($rAnnotation->{"2L"}[1]{feat}, qr/1/, "ncRNA, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[1]{feat}, qr/a/, "pre_miRNA, 2L");

like($rAnnotation->{"2L"}[2]{feat}, qr/e/, "exon, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[2]{feat}, qr/T/, "three_prime_UTR, 2L");

like($rAnnotation->{"2L"}[3]{feat}, qr/i/, "intron, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/C/, "CDS, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/T/, "three_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/h/, "three_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/F/, "five_prime_UTR, 2L");
like($rAnnotation->{"2L"}[3]{feat}, qr/p/, "pseudogene, 2L");

like($rAnnotation->{"2R"}[1]{feat}, qr/4/, "snRNA, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[1]{feat}, qr/a/, "pre_miRNA, 2R");

like($rAnnotation->{"2R"}[2]{feat}, qr/2/, "tRNA, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[2]{feat}, qr/t/, "transposable_elenemt, 2R");

like($rAnnotation->{"2R"}[3]{feat}, qr/i/, "intron, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/T/, "three_prime_UTR, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[3]{feat}, qr/t/, "transposable_elenemt, 2R");

like($rAnnotation->{"2R"}[4]{feat}, qr/5/, "rRNA, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/F/, "five_prime_UTR, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[4]{feat}, qr/a/, "pre_miRNA, 2R");

like($rAnnotation->{"2R"}[5]{feat}, qr/3/, "snoRNA, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/m/, "miRNA, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/h/, "enhancer, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/r/, "regulatory_region, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/p/, "pseudogene, 2R");
like($rAnnotation->{"2R"}[5]{feat}, qr/a/, "pre_miRNA, 2R");

#
# load_GTF_file_CDS_with_gene_IDs
#

can_ok('GffGtfParser', 'load_GTF_file_CDS_with_gene_IDs');

$gtf = "chr4\tdm3_flyBaseGene\texon\t264260\t264374\t0.000000\t+\t.\tgene_id \"CG1674-RB\"; transcript_id \"CG1674-RB\";\n".
	   "chr4\tdm3_flyBaseGene\tCDS\t2\t4\t0.000000\t\+\t2\tgene_id \"CG1674-RB\"; transcript_id \"CG1674-RB\";"; 


$rAnnotation={};
$ptrGeneIDs={};
$rAnnotation = load_GTF_file_CDS_with_gene_IDs(\$gtf, $ptrGeneIDs);
is($rAnnotation->{4}[2]{feat},"C", "load only CDS");
is($rAnnotation->{4}[3]{feat},"C", "load only CDS");
is($rAnnotation->{4}[4]{feat},"C", "load only CDS");
is(scalar keys %$ptrGeneIDs, 1, "load only CDS, geneIDs size ok");
#print Dumper($ptrGeneIDs);
is($ptrGeneIDs->{"CG1674-RB"}[0]{chromosome},4,"load only CDS, geneID content ok");
is($ptrGeneIDs->{"CG1674-RB"}[0]{start},2,"load only CDS, geneID content ok");
is($ptrGeneIDs->{"CG1674-RB"}[0]{end},4,"load only CDS, geneID content ok");

#
# get_characteristics_of_genome_gff_pileup 
#

can_ok('GffGtfParser', 'get_characteristics_of_genome_gff_pileup');

$rNumbers={};

$rNumbers = get_characteristics_of_genome_gff_pileup(
\$gffAllFeatures, \$PILEUP_FILE,
$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
$POOL_SIZE);

is($rNumbers->{CDS}{pi}, 0.750105011201128,"pi");
is($rNumbers->{CDS}{theta}, 0.689766164273035,"theta");
is($rNumbers->{CDS}{D}, 0.684940610557322,"D"); 
