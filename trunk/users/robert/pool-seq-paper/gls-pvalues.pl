#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path;
# do enter correct path here
use MaxCoverage;
use Synchronized;
use SynchronizeUtility;
use MajorAlleles; # get the two major allele
use Test;

# Author: Robert Kofler
# Define the variables
my $input;
my $output="";
my $help=0;
my $test=0;
my $verbose=1;

my $mincount=2;
my $mincoverage=4;
my $usermaxcoverage;
my $minlogpvalue=0.0;
my $removetemp=0;

# --input /Users/robertkofler/pub/PoPoolation2/Walkthrough/demo-data/cmh/small-test.sync --output /Users/robertkofler/pub/PoPoolation2/Walkthrough/demo-data/cmh/small-test.cmh --population 1,2,3,4 --min-count 2 --min-coverage 4 --max-coverage 200

GetOptions(
    "input=s"	    =>\$input,
    "output=s"	    =>\$output,
    "min-count=s"   =>\$mincount,
    "min-coverage=i"=>\$mincoverage,
    "max-coverage=s"=>\$usermaxcoverage,
    "test"          =>\$test,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
die "Tests currently not supported" if $test;
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $output;
#pod2usage(-msg=>"Minimum coverage must be equal or larger than minimum count",-verbose=>1) unless $mincoverage>= $mincount;
pod2usage(-msg=>"Maximum coverage has to be provided",-verbose=>1) unless $usermaxcoverage;




my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $syncparser=get_sumsnp_synparser($mincount,$mincoverage,$maxcoverage);

my $rinput=$output.".rin";
my $routput=$output.".rout";

print "Reading sync file and writing temporary R output file\n";
CMHUtil::write_Rinput($input,$rinput,$syncparser);

print "Calling R, to calculate the gls test statistic\n";
system("R --vanilla --slave <$rinput >$routput");

print "Parsing R-output and writing output file\n";
CMHUtil::write_output($routput,$output);

if($removetemp)
{
	print "Removing temporary files\n";
	unlink($rinput);
	unlink($routput);
}
print "Done\n";

exit(0);



{
	package CMHUtil;
	use strict;
	use warnings;
	use List::Util qw[min max];
	use MaxCoverage;
	use Synchronized;
	use SynchronizeUtility;
	use MajorAlleles; # get the two major allele
	
	sub write_output
	{
		my $routput=shift;
		my $output=shift;

		
		open my $ifh,"<", $routput or die "Could not open input file\n";
		open my $ofh,">",$output or die "Could not open output file\n";
		
		while(1)
		{
			#[1] "2R\t2296\tN\t90:10:0:0:0:0\t100:0:0:0:0:0\t100:0:0:0:0:0\t100:0:0:0:0:0"
			#[1] 0.003583457
			my $line=<$ifh>;
			last unless $line;
			my $dispersion=<$ifh>;
			my $pvalue_lib=<$ifh>;
			my $pvalue_rep=<$ifh>;
			chomp $line; chomp $dispersion; chomp $pvalue_lib; chomp $pvalue_rep;
			$line=~s/^\S+\s//;
			$line=~s/^"//;
			$line=~s/"$//;
			$line=~s/\\t/\t/g;
			$dispersion=~s/^\S+\s//;
			$pvalue_lib=~s/^\S+\s//;
			$pvalue_rep=~s/^\S+\s//;
			#$pvalue="1.0" if $pvalue eq "NaN"; 	# stupid mantelhaenszeltest prodcues NaN for example mantelhaen.test(array(c(100,100,0,0,100,100,0,0,100,100,0,0),dim=c(2,2,3)),alternative=c("two.sided"))
								# this is clearly no differentiation thus 1.0 (necessary as it fucks up sorting by significance)

			print $ofh $line."\t".$dispersion."\t".$pvalue_lib."\t".$pvalue_rep."\n";
		}
		close $ofh;
		close $ifh;
	}
	
	

	
	sub _write_common
	{
		my $fh=shift;
print $fh <<PERLSUCKS;
getdispersion<-function(counts)
{
  lib =rep(c("a","b","c","d"),each=4)
  rep =rep(c("r1","r1","r2","r2"),4)
  allele =rep(c("a1","a2"),8)
  data<-data.frame(lib,rep,allele,counts)
  
  model1a= glm(counts~lib*rep+allele*lib,data=data,family=quasipoisson)
  model1b= glm(counts~lib*rep+allele*rep,data=data,family=quasipoisson)
  model3= glm(counts~lib*rep+allele*lib + allele*rep,data=data,family=quasipoisson)

  a1<-anova(model1a,model3,test="LRT") # replicate Influence
  a2<-anova(model1b,model3,test="LRT") # library Influence
  
  dispersion<-summary(model3)\$dispersion
  p.replicate <- a1\$Pr[2]
  p.library <- a2\$Pr[2]
  return(c(dispersion,p.library,p.replicate))
}
PERLSUCKS

	}
	
	sub write_Rinput
	{
		my $syncfile=shift;
		my $rinput=shift;
		my $syncparser=shift;
		
		open my $ifh, "<", $syncfile or die "Could not open input file";
		open my $ofh, ">", $rinput or die "Could not open routput file";
		_write_common($ofh);
		while(my $line=<$ifh>)
		{
			chomp $line;
			my $e=$syncparser->($line);
			next unless $e->{ispuresnp};
			_write_single_entry($ofh,$line,$e);
		}
		close $ofh;
		close $ifh;
	}
	
	
	

	
	
	
	
	

	
	
	sub _write_single_entry
	{


#lib rep allele counts exp.counter
#1    a  r1     a1     36           a
#2    a  r1     a2     26           a
#3    a  r2     a1     25           b
#4    a  r2     a2     26           b
#5    b  r1     a1     31           c
#6    b  r1     a2     34           c
#7    b  r2     a1     29           d
#8    b  r2     a2     27           d
#9    c  r1     a1     34           e
#10   c  r1     a2     26           e
#11   c  r2     a1     25           f
#12   c  r2     a2     41           f
#13   d  r1     a1     35           g
#14   d  r1     a2     29           g
#15   d  r2     a1     31           h
#16   d  r2     a2     25           h
# l1-r1-a1, l1-r1-a2, l1-r2-a1,l1-r2-a2....

		my $ofh=shift;
		my $line=shift;
		my $entry=shift;


			my $countstring = _get_countstring($entry);
			print $ofh <<PERLSUCKS;
print("$line")
tmp<-getdispersion($countstring)
print(tmp[1])
print(tmp[2])
print(tmp[3])
PERLSUCKS


	}
	
	sub _get_countstring
	{
		my $e=shift;
		my ($major,$minor) = MajorAlleles::get_major_minor_alleles($e->{samples});
		my $as=[];
		for my $s (@{$e->{samples}})
		{
			push @$as,$s->{$major};
			push @$as,$s->{$minor};
		}
		my $cs=join(",",@$as);
		my $countstr="c($cs)";
		return $countstr;
	}
	
	
	

	
	
	
    
}








=head1 NAME

gls-dispersion.pl

=head1 SYNOPSIS

 perl gls-dispersion.pl --input ip.sync --output dispersion.txt --min-count 1 --min-coverage 10 --max-coverage 50

=head1 AUTHORS

 Robert Kofler
 Christian Schloetterer

=cut
