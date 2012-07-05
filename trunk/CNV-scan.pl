#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use List::Util qw[min max];



# --pool-size 500 --min-count 2 --min-coverage 4 --window-size 1000 --step-size 1000 --input test/snp.merge --output test/test.fst


my $input;
my $output="";
my $help=0;
my $test=0;
my $mincoverage=4;
my $compare="";
my $threshold=1.5;
my $minLength=100;

GetOptions(
    "input=s"	        =>\$input,
    "output=s"          =>\$output,
    "compare=s"         =>\$compare,
    "threshold=f"       =>\$threshold,
    "min-coverage=i"    =>\$mincoverage,
    "min-length=i"      =>\$minLength,
    "test"              =>\$test,
    "help"	        =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

# too many arguments should result in an error
pod2usage(-msg=>"Wrong options",-verbose=>1) if @ARGV;
pod2usage(-verbose=>2) if $help;
Test::runTests() if $test;

pod2usage(-msg=>"Input file does not exist",-verbose=>1) unless -e $input;
pod2usage(-msg=>"No output file has been provided",-verbose=>1) unless $output;
pod2usage(-msg=>"Please choose which populations to compare",-verbose=>1) unless $compare;
pod2usage(-msg=>"Minimum coverage <1 not allowed",-verbose=>1) if $mincoverage<1;
pod2usage(-msg=>"Threshold has to be larger than zero",-verbose=>1) unless $threshold>1;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using compare\t$compare\n";
print $pfh "Using min-coverage\t$mincoverage\n";
print $pfh "Using min-length\t$minLength\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my($pop1,$pop2) = split /:/,$compare;
my($avcov1,$avcov2)=Utility::get_average_coverage($input,$pop1,$pop2);

open my $ifh, "<", $input or die "Could not open input file";
open my $ofh, ">", $output or die "Could not open output file";

my $hcnv=Utility::get_hsp_finder($mincoverage,$threshold,$avcov1,$avcov2);
my $lcnv=Utility::get_hsp_finder($mincoverage,$threshold,$avcov2,$avcov1);

print "Start computing the CNV's between population $pop1 and population $pop2\n";
while(my $line = <$ifh>)
{
    chomp $line;
    my $e=Utility::parseSync($line);
    my($chr,$pos)=($e->{chr},$e->{pos});
    my ($c1,$c2)=($e->{data}[$pop1-1]{totcov},$e->{data}[$pop2-1]{totcov});
    my $hhs=$hcnv->($chr,$pos,$c1,$c2);
    my $lhs=$lcnv->($chr,$pos,$c2,$c1);
    
    Utility::print_cnv($ofh,$minLength,$hhs,$lhs);
}

my $hhs=$hcnv->("",0,0,0,1);
my $lhs=$lcnv->("",0,0,0,1);
Utility::print_cnv($ofh,$minLength,$hhs,$lhs);

print "Finished..\n";

exit;


{
    package Utility;
    use strict;
    use warnings;
    
    
    sub print_cnv
    {
        my $fh=shift;
        my $minLength=shift;
        my $hhs=shift;
        my $lhs=shift;
        
        if($hhs)
        {
            if($hhs->{cov_count}>=$minLength)
            {
                # chr, start, end, length, score, avcnv
                my $covf=sprintf("%.3f",$hhs->{cov_frac});
                my $avcnv=sprintf("%.3f",$hhs->{avcnv});
                my $score=sprintf("%.3f",$hhs->{score});
                print $fh "$hhs->{chr}\t$hhs->{start}\t$hhs->{end}\t$hhs->{length}\t$covf\t$avcnv\t$score\n";
            }
        }
        if($lhs)
        {
            if($lhs->{cov_count}>=$minLength)
            {
                my $covf=sprintf("%.3f",$lhs->{cov_frac});
                my $avcnv=sprintf("%.3f",$lhs->{avcnv});
                my $score=sprintf("%.3f",$lhs->{score});
                print $fh "$lhs->{chr}\t$lhs->{start}\t$lhs->{end}\t$lhs->{length}\t$covf\t-$avcnv\t$score\n";
            }
        }
        
    }
    
    
    sub get_hsp_finder
    {
        my $mincov=shift;
        my $threshold=shift;
        my $avc1=shift;
        my $avc2=shift;
        
        my $t={};
        
        $t=_reset_t($t);
        

        return sub
        {

            my $chr=shift;
            my $pos=shift;
            my $c1=shift;
            my $c2=shift;
            my $retflag= shift or 0;
            return _return_hs($t,$threshold) if $retflag;
            
            $t->{prevchr}=$chr unless $t->{prevchr};
            $t->{prevpos}=$pos-1 unless $t->{prevpos};
            

            # define a high score;
            my $hs=undef;
            if($chr ne $t->{prevchr})
            {
                $hs=_return_hs($t,$threshold);
                $t=_reset_t($t);
            }
            
            # set the pos score to zero if the coverage is not sufficient
            # if coverage is sufficient increase the counter
            my $pos_score=0;
            if($c1>=$mincov and $c2>=$mincov)
            {
                my $norm_c1=$c1/$avc1;
                my $norm_c2=$c2/$avc2;
                my $ratio=$norm_c1/$norm_c2;
                $pos_score=$ratio-$threshold;
                $t->{run_cov_count}++;
            }

            $t->{run_score}+=$pos_score;

            
            if($t->{run_score} > $t->{hs_score})
            {
                $t->{hs_start}      =$pos if($t->{hs_start}==-1);
                $t->{hs_score}      =$t->{run_score};
                $t->{hs_cov_count}  =$t->{run_cov_count};
                $t->{hs_end}        =$pos;
            }
            elsif($t->{run_score}<0)
            {
                $hs=_return_hs($t,$threshold);
                $t=_reset_t($t);
            }
            return $hs;
        }
    }
    
    sub _reset_t
    {
        my $t=shift;
        
        $t={
          prevchr=>"",
          prevpos=>"",
          hs_start=>-1,
          hs_end=>-1,
          hs_score=>0,
          hs_cov_count=>0,
          run_score=>0,
          run_cov_count=>0
        };
        return $t;
        
    }
    
    
    sub _return_hs
    {
        my $t=shift;
        my $threshold=shift;
        
        return undef unless $t->{hs_score}>0;
        
        
        # chr, start, end, length, score, avcnv
        my $hs=
        {
            start=>$t->{hs_start},
            end=>$t->{hs_end},
            length=>$t->{hs_end}-$t->{hs_start}+1,
            cov_count=>$t->{hs_cov_count},
            chr=>$t->{prevchr},
            score=>$t->{hs_score}
        };
        
        $hs->{avcnv}    =($hs->{score}/$t->{hs_cov_count})+$threshold;
        $hs->{cov_frac} =$hs->{cov_count}/$hs->{length};
        return $hs;
        
    }
    
    
    
    sub get_average_coverage
    {
        my $input=shift;
        my $pop1=shift;
        my $pop2=shift;
        
        open my $ifh,"<",$input or die "could not open file handle";
        
        my $p1c=0;
        my $p2c=0;
        my $count=0;
        
        print "Computing average coverages for population $pop1 and population $pop2\n";
        print "Start parsing the synchronized file\n";
        
        while(my $l=<$ifh>)
        {
            chomp $l;
            my $e=parseSync($l);
            my($c1,$c2)=($e->{data}[$pop1-1]{totcov},$e->{data}[$pop2-1]{totcov});
            $p1c+=$c1;
            $p2c+=$c2;
            $count++;
        }
        
        my $acp1=$p1c/$count;
        my $acp2=$p2c/$count;
        print "Finished computing the average coverage\n";
        print "base pairs covered: $count\n";
        print "Average coverage for population $pop1: $acp1\n";
        print "Average coverage for population $pop2: $acp2\n";
        return ($acp1,$acp2); 
    }
    
    
    
    
    

    
    sub parseSync
    {
        my $line=shift;
        
        #gi|42519920|ref|NC_002978.6     1       T       0:15:0:0:0:0    0:34:0:0:0:0    0:2:0:0:0:0     0:7:0:0:0:0     0:9:0:0:0:0
        #gi|42519920|ref|NC_002978.6     2       G       0:0:0:15:0:0    0:0:0:34:0:0    0:0:0:2:0:0     0:0:0:7:0:0     0:0:0:9:0:0

        my @a=split /\t/,$line;
        my $chr=shift @a;
        my $pos=shift @a;
        my $rc=shift @a;
        
        my @data=();
        foreach my $e (@a)
        {
            my $p=_parse_entry($e);
            push @data,$p;
        }
        
        return
        {
            chr=>$chr,
            pos=>$pos,
            rc=>$rc,
            data=>\@data
        };
    }


    sub _parse_entry
    {
        my $entry  = shift;
        
        my $toret;
        if($entry=~/-/)
        {
            $toret={
                eucov=>0,
                totcov=>0,
                A=>0,
                T=>0,
                C=>0,
                G=>0,
                N=>0,
                del=>0
                };
            }
            else
            {
                my @parts=split /:/,$entry;
                
                die "failed parsing $entry; does not have the correct number of entries" unless @parts ==6;
                
                $toret = {
                     A  => $parts[0],
                     T  => $parts[1],
                     C  => $parts[2],
                     G  => $parts[3],
                     N  => $parts[4],
                     del=> $parts[5]
                    }; 
                $toret->{eucov}  = ($toret->{A} + $toret->{T} + $toret->{C} + $toret->{G});
                $toret->{totcov} = ($toret->{eucov}+ $toret->{N} + $toret->{del} );
                # A, T, C, G, N, del, eucov, totcov, valid_cov (totcov), a_desc
            }
    return $toret;
   }
}



=head1 NAME

CNV-scan.pl - Calculate copy number variations (CNV) between two pooled popultions

=head1 SYNOPSIS

 CNV-scan.pl --input populations.sync --output cnv.txt --compare 2:3 --threshold 1.5 --min-coverage 4 --min-length 100

=head1 OPTIONS

=over 4

=item B<--input>

The input file. Has to be synchronized pileup file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--compare>

The samples which should be compared. A synchronized file may contain >=2 populations. The script can only compare two populations at a time. example: "--compare 2:3"; mandatory parameter

=item B<--threshold>

CNVs to report should on average have a fold difference of <--threshold>; default=1.5

=item B<--min-coverage>

the minimum coverage; used for CNV identification, the coverage in ALL populations has to be higher or equal to this threshold, otherwise no SNP will be called. default=4

=item B<--min-length>

the minimum length of CNV; only regions meeting the <--min-coverage> will be considered. default=100

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS

Combine different pileup files

=head2 INPUT

A synchronized file, for example

 Unknown_group_104	5943	N	0:0:10:0:0:0	0:0:10:0:0:0	0:0:10:0:0:0
 Unknown_group_104	5944	N	0:8:0:0:0:0	0:8:0:0:0:0	0:8:0:0:0:0
 Unknown_group_104	5945	N	0:0:0:8:0:0	0:0:0:8:0:0	0:0:0:8:0:0
 Unknown_group_104	5946	N	0:0:9:0:0:0	0:0:9:0:0:0	0:0:9:0:0:0
 Unknown_group_104	5947	N	0:7:0:0:0:0	0:7:0:0:0:0	0:7:0:0:0:0
 
 col1: reference contig (chromosome)
 col2: position in the reference contig
 col3: reference character
 col4: population 1
 col5: population 2
 coln: population n
 
 population data are in the form
 A:T:C:G:N:*
 A: count of character A
 T: count of character T
 C: count of character C
 G: count of character G
 N: count of character N
 *: deletion, count of deletion

Only the characters A,T,C,G are considered for the coverage.
The SNP site is ignored if a deletion is found in any population

=head2 OUTPUT <output>

 2R     86259   87202   944     1.000   -1.770  254.672
 2R     89089   90298   1210    0.840   2.877   1400.364

 col1: chromosome of the CNV
 col2: start position of the CNV
 col3: end position of the CNV
 col4: length of the CNV
 col5: fraction of the CNV which is sufficiently covered; positions with insufficient coverage are ignored!!
 col6: average fold copy number difference between the two populations; positive if the first population has more copies ( first with respect to <--compare>) negative if the second population has more copies
 col7: score of the CNV

=head1 AUTHORS

Robert Kofler

Christian Schloetterer

=cut



