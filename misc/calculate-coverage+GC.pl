#!/usr/bin/perl-w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use DBI;
use DBI::DBD;
# Program name: calculate-coverage+GC.pl
# Author: Ram Vinay Pandey 
# Date: 09-09-2010
# Description: This script calculates average coverage and GC,AT content for sliding window as well as for whole chromosome.
   
# Define the variables
my $input;
my $output;
my $windowsize=10000;
my $step=1000;
my $minCoverageFraction=0.6;

my $help=0;
my $test=0;
my $verbose=1;


my $usage="perl $0 --input test.pileup --output coverage+GC-output.txt --window-size 10000 --step-size 1000 --min-covered-fraction 0.6 \n";

GetOptions(
    "input=s"	    =>\$input,
    "output=s"	    =>\$output,
    "window-size=s"  =>\$windowsize,
    "step-size=i"  =>\$step,
    "min-covered-fraction=f"=>\$minCoverageFraction,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
Test::runTests() if $test;

  
pod2usage(-msg=>"\n\t$usage\n\tA pileup input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"\n\t$usage\n\tProvide output file name\n",-verbose=>1) unless $output;


print "Chromosome wise Started at: \t". localtime() .".....\n";

my ( $name, $path, $extension ) = File::Basename::fileparse ( $output, '\..*' );

my ($window_average_cov_file,$chromosome_average_cov_file) = ("","");

$window_average_cov_file = $path.$name."-window-average-coverage-GC-content.txt";
$chromosome_average_cov_file = $path.$name."-chromosome-average-coverage-GC-content.txt";


open my $ofh2, ">$chromosome_average_cov_file" or die "Could not open output file";

Coverage::_caculate_Chr_Coverage_GC($input,$ofh2);

print "Chromosome wise Completed at: \t". localtime() .".....\n\n";

print "Sliding window wise Started at: \t". localtime() .".....\n";
open my $ofh1, ">$window_average_cov_file" or die "Could not open output file";

#=cut
my $reader;

$reader=BpSlider->new($input,$windowsize,$step);


while(my $window=$reader->nextWindow())
{
        my $chr=$window->{chr};
	my $start=$window->{start};
	my $end=$window->{end};
        my $win=$window->{window};
	my $pos=$window->{middle};
	my $above=$window->{count_covered};
        my $data=$window->{data};
	
        next unless @$data;
	#print "$chr\t$pos\n";
	Coverage::_caculateCoverage_GC($data,$pos,$chr,$win,$minCoverageFraction,$ofh1);
	#Coverage::_calculateGC($data,$pos,$chr,$win,$minCoverageFraction,$ofh2);

}

print "Sliding window wise Completed at: \t". localtime() .".....\n";



exit;




{
    use warnings;
    use strict;
    package BpSlider;

    sub new
    {
        my $class=shift;
        my $file=shift;
        my $window=shift;
        my $step=shift;
        
        open my $fh,"<$file" or die "Could not open file handle";
        
        return bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            file=>$file,
            fh=>$fh,
            curwin=>[],
            buffer=>[]
        },__PACKAGE__;
    }
    
    sub nextWindow
    {
        my $self=shift;
	my $file=$self->{file}; 
	
        #get the current window, and the current chromosome
        my $curwin=$self->{curwin};
        
        my $curChr="";
        $curChr=$curwin->[0]{chr} if @$curwin;
        
        my $resetchr=0;
        
        # empty unnecessary entries
        EMPTY: while(@$curwin)
        {
            my $e=shift @$curwin;
            if($e->{pos}>$self->{lower})
            {
                unshift @$curwin, $e;
                last EMPTY;
            }
            
        }
        
        # fill with novel entries
        my $line;
        FILL:while($line=$self->_nextline)
        {
            my $e=Utility::_parseLightwight($line);
            $curChr=$e->{chr} unless $curChr;
            
            
            if($e->{chr} eq $curChr && $e->{pos} <= $self->{upper})
            {
                push @$curwin,$e;
            }
            else
            {
                $resetchr=1 if $e->{chr} ne $curChr;
                $self->_bufferline($line);
                last FILL;
            }
        }
        
        return undef unless $curChr;
        
        
        my $toret=Utility::_annotateWindow($curwin,$curChr,$self->{lower},$self->{upper},$self->{window});
        
        if($resetchr or not defined($line))
        {
            # we transgressed the boundaries to the next chromosome
            # reset the windows and the current buffer
            $self->{lower}=0;
            $self->{upper}=$self->{window};
            $self->{curwin}=[];
        }
        else
        {
            # next time we will still be in the same chromosome
            # increase the upper and lower boundaries by the stepsize and set the current buffer
            $self->{upper}+=$self->{step};
            $self->{lower}+=$self->{step};
            $self->{curwin}=$curwin;
        }

        return $toret;
    }
    
    
    
    sub _nextline
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub _bufferline
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }
    
    
}



{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    
    sub _annotateWindow
    {
        my $curwin=shift;
        my $chr=shift;
        my $start=shift;
        my $end=shift;
        my $window=shift;

        my $snps=0;
        my $aboveCoverage=0;
        foreach(@$curwin)
        {
            $aboveCoverage++ if $_->{iscov};
        }

        return
        {
            chr=>$chr,
            start=>$start,
            end=>$end,
	    middle=>int(($end+1+$start)/2),
            count_covered=>$aboveCoverage,
            window=>$window,
            data=>$curwin      
        };
    }
    
    
    
    sub _parseLightwight
    {
        my $line=shift;
	
        chomp $line;
        my @a=split /\s+/,$line;
        my $chr=shift @a;
        my $pos=shift @a;
        my $rc=shift @a;
        my $cov=shift @a;
	my $seq=shift @a;
	my $qual=shift @a;
	        
        
        my $en={
            chr=>$chr,
            pos=>$pos,
            refchar=>$rc,
            cov=>$cov,
            seq=>$seq,
	    qual=>$qual,
	    iscov=>1
        };        
        return $en;
    }

    
    sub _getChromosome {
	my $file=shift;
	my %chromosome = ();
	my $ct=0;
	
	open (PILEUP,"<$file") || die ("could not open pileup file for reading $!\n");
	while (<PILEUP>)
	{
	    chomp;
	    my $line = $_;
	    # discard blank line
	    if ($line =~ m/^\s*$/g) {
		next;
	    }
	    else {
		$ct++;
		my @a=split /\s+/,$line;
		my $chr = shift @a;
		$chromosome{$chr}=$ct;
	    }
	}
	
	close PILEUP;
	return %chromosome;
    }
    
    ##########################################################################################
    # return_value function return all columns for each line (\t seperated) of input file. 
    ##########################################################################################
    
    sub _returnValue {
	    my ($fh,$priv_pos,$priv_chr,$priv_rc,$priv_cov) = (shift,shift,shift,shift,shift);
	    
	    # Declaration and initialization of variables
	    my $curline = "";
	    my @vals = ();
	    my %vals = ();
	    
	    
	    if (defined ($curline = <$fh>)) {
	    	$curline = $curline;
	    }
	    else {
		    $curline = "";
	    }

	    #2L_mau_20N      45      T       33      .................................       bbba_babWaaa\ba`abbbbbaa`bbabbaaa
	    my @colnames = ("chr", "pos", "refnuc", "cov", "nucsori", "qual");
	    
	    if ($curline ne "") {
		chomp $curline;
		@vals = split("\t",$curline);
		@vals{@colnames} = @vals;

		$vals{"privchr"} = $priv_chr;
	    	$vals{"privpos"} = $priv_pos;
		$vals{"privrefnuc"} = $priv_rc;
	    	$vals{"privcov"} = $priv_cov;
		
		$vals{"diff"} = abs($vals{"pos"} - $priv_pos);
	    	$vals{"fh"} = $fh;
	    	$vals{"line"} = $curline;
	    }
	    else {
		@vals = ("", 0, "", "", "", "");
		@vals{@colnames} = @vals;

		$vals{"privchr"} = $priv_chr;
		$vals{"privpos"} = $priv_pos;
		$vals{"privrefnuc"} = $priv_rc;
	    	$vals{"privcov"} = $priv_cov;
		$vals{"diff"} = abs($vals{"pos"} - $priv_pos);
	    	$vals{"fh"} = $fh;
	    	$vals{"line"} = $curline;
	    }

	    return \%vals
    
    }
    
    
}




{
	package Coverage;
	use strict;
	use warnings;
	
	sub _caculateCoverage_GC {
		
		my $data=shift;
		my $pos=shift;
		my $chr=shift;
		my $win=shift;
		my $minCovFraction=shift;
		my $ofh=shift;
		
		my $covered_fraction=0;
		my $considered_pos = 0;
		
		my ($total_cov,$average_cov) = (0,0);
		my ($act,$tct,$gct,$cct,$nct) = (0,0,0,0,0);
		my ($gc_total,$at_total,$gc_percent,$at_percent) = (0,0,0,0);
		
		$considered_pos = @$data;
		if ($considered_pos>0) {
			$covered_fraction = $considered_pos/$win;
			$covered_fraction = sprintf "%.2f",$covered_fraction;
		}
		
		foreach my $d (@$data) {
			
			$total_cov = $total_cov+$d->{cov};
			
			if ($d->{refchar} =~ /A/i) {
				$act++;
			}
			elsif ($d->{refchar} =~ /T/i) {
				$tct++;
			}
			elsif ($d->{refchar} =~ /G/i) {
				$gct++;
			}
			elsif ($d->{refchar} =~ /C/i) {
				$cct++;
			}
			elsif ($d->{refchar} =~ /N/i) {
				$nct++;
			}

		}
		
		
		
		$gc_total = $gct+$cct;
		$at_total = $act+$tct;
		
		
		if(($considered_pos>0) and ($gc_total>0)) {
			$gc_percent = ($gc_total/$considered_pos)*100;
			#$gc_percent = sprintf "%.2f",$gc_percent;
		}
		else {
			$gc_percent = "0.0000";
		}
		
		if(($considered_pos>0) and ($at_total>0)) {
			$at_percent = ($at_total/$considered_pos)*100;
			$at_percent = sprintf "%.4f",$at_percent;
		}
		else {
			$at_percent = "0.0000";
		}
		
		if(($considered_pos>0) and ($total_cov>0)) {
			$average_cov = ($total_cov/$considered_pos);
			$average_cov = sprintf "%.4f",$average_cov;
		}
		else {
			$average_cov = "0.0000";
		}
		
		if ($covered_fraction>=$minCovFraction) {
			print $ofh "$chr\t$pos\t$covered_fraction\t$average_cov\t$gc_percent\t$at_percent\n";
		}
		else {
			print $ofh "$chr\t$pos\t$covered_fraction\tna\tna\tna\n";
		}
		
		
		
	}
	
	sub _caculate_Chr_Coverage_GC {
		my $file=shift;
		my $ofh=shift;
		
		my ($priv_pos,$ct,$lineval,$priv_chr,$priv_rc,$priv_cov) = (0,0,"","",0);
		
		my ($total_cov,$average_cov) = (0,0);
		my ($act,$tct,$gct,$cct,$nct) = (0,0,0,0,0);
		my ($gc_total,$at_total,$gc_percent,$at_percent) = (0,0,0,0);
		my $considered_pos=0;
		
		open my $fh,"<$file" or die "Could not open pileup input file $file for reading $!\n";
		
		while(1) {
	    
			$lineval = Utility::_returnValue($fh,$priv_pos,$priv_chr,$priv_rc,$priv_cov);
			
			$priv_chr = $lineval->{chr};
			$priv_pos = $lineval->{pos};
			$priv_rc = $lineval->{refnuc};
			$priv_cov = $lineval->{cov};
			$fh = $lineval->{fh};
			if (($lineval->{pos}>0) and (($lineval->{privpos}>0))) {
			    
			    if ($lineval->{privchr} eq $lineval->{chr}) {
				
				$considered_pos++;
				$total_cov = $total_cov+$lineval->{privcov};
				
				if ($lineval->{privrefnuc} =~ /A/i) {
				$act++;
				}
				elsif ($lineval->{privrefnuc} =~ /T/i) {
					$tct++;
				}
				elsif ($lineval->{privrefnuc} =~ /G/i) {
					$gct++;
				}
				elsif ($lineval->{privrefnuc} =~ /C/i) {
					$cct++;
				}
				elsif ($lineval->{privrefnuc} =~ /N/i) {
					$nct++;
				}
				
				#print "0: $lineval->{privchr}\t$lineval->{privpos}\t$lineval->{privcov}\t$lineval->{privrefnuc}\t$lineval->{chr}\t$lineval->{pos}\t$lineval->{cov}\t$lineval->{refnuc}\n";

				
			    }

			}

			
			## Write average coverage GC and AT content
			if ($lineval->{privpos}>0) {
			    
			    if ($lineval->{privchr} ne $lineval->{chr}) {
		
				$considered_pos++;
				#print "2: $lineval->{privchr}\t$lineval->{privpos}\t$lineval->{privcov}\t$lineval->{privrefnuc}\t$lineval->{chr}\t$lineval->{pos}\t$lineval->{cov}\t$lineval->{refnuc}\n";
				
				$total_cov = $total_cov+$lineval->{privcov};
				
				if ($lineval->{privrefnuc} =~ /A/i) {
				$act++;
				}
				elsif ($lineval->{privrefnuc} =~ /T/i) {
					$tct++;
				}
				elsif ($lineval->{privrefnuc} =~ /G/i) {
					$gct++;
				}
				elsif ($lineval->{privrefnuc} =~ /C/i) {
					$cct++;
				}
				elsif ($lineval->{privrefnuc} =~ /N/i) {
					$nct++;
				}
				
				
				$gc_total = $gct+$cct;
				$at_total = $act+$tct;
				
				
				if(($considered_pos>0) and ($gc_total>0)) {
					$gc_percent = ($gc_total/$considered_pos)*100;
					#$gc_percent = sprintf "%.2f",$gc_percent;
				}
				else {
					$gc_percent = "0.00";
				}
				
				if(($considered_pos>0) and ($at_total>0)) {
					$at_percent = ($at_total/$considered_pos)*100;
					#$at_percent = sprintf "%.2f",$at_percent;
				}
				else {
					$at_percent = "0.00";
				}
				
				if(($considered_pos>0) and ($total_cov>0)) {
					$average_cov = ($total_cov/$considered_pos);
					#$average_cov = sprintf "%.2f",$average_cov;
				}
				else {
					$average_cov = "0.00";
				}
				
					print $ofh "$lineval->{privchr}\t$average_cov\t$gc_percent\t$at_percent\n";

				### empty all variables
				($total_cov,$average_cov) = (0,0);
				($act,$tct,$gct,$cct,$nct) = (0,0,0,0,0);
				($gc_total,$at_total,$gc_percent,$at_percent) = (0,0,0,0);
				$considered_pos=0;
				#$ct1=0;
			    }
			}
		
			#my $diff = 
			last if ($lineval->{line} eq "");
		    }

	}
	sub _calculateGC {
		
		my $data=shift;
		my $pos=shift;
		my $chr=shift;
		my $win=shift;
		my $minCovFraction=shift;
		my $ofh=shift;
		
		my $covered_fraction=0;
		my $considered_pos = 0;
		
		$considered_pos = @$data;
		if ($considered_pos>0) {
			$covered_fraction = $considered_pos/$win;
			$covered_fraction = sprintf "%.2f",$covered_fraction;
		}
		
		my ($act,$tct,$gct,$cct,$nct) = (0,0,0,0,0);
		
		foreach my $d (@$data) {
			
			if ($d->{refchar} =~ /A/i) {
				$act++;
			}
			elsif ($d->{refchar} =~ /T/i) {
				$tct++;
			}
			elsif ($d->{refchar} =~ /G/i) {
				$gct++;
			}
			elsif ($d->{refchar} =~ /C/i) {
				$cct++;
			}
			elsif ($d->{refchar} =~ /N/i) {
				$nct++;
			}
		}
		
		my ($gc_total,$at_total,$gc_percent,$at_percent) = (0,0,0,0);
		$gc_total = $gct+$cct;
		$at_total = $act+$tct;
		
		
		if(($considered_pos>0) and ($gc_total>0)) {
			$gc_percent = ($gc_total/$considered_pos)*100;
			$gc_percent = sprintf "%.2f",$gc_percent;
		}
		else {
			$gc_percent = "0.00";
		}
		
		if(($considered_pos>0) and ($at_total>0)) {
			$at_percent = ($at_total/$considered_pos)*100;
			$at_percent = sprintf "%.2f",$at_percent;
		}
		else {
			$at_percent = "0.00";
		}
		
		if ($covered_fraction>=$minCovFraction) {
			
		}
		else {
			
		}
	
		
	}
	
}



=head1 NAME

calculate-coverage+GC.pl - This script takes a pileup file as input and calculates average coverage, GC and AT content 1) for sliding windows along each chromosome and 2) across each entire chromosome.

=head1 SYNOPSIS

 calculate-coverage+GC.pl --input dmel-pileup.pileup --output coverage.txt --window-size 10000 --step-size 1000 --min-covered-fraction 0.6 

=head1 OPTIONS

=over 4

=item B<--input>

The input file in pileup format. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--min-covered-fraction>

the minimum fraction of a window being covered by at least one read; float; default=0.6


=item B<--window-size>

the size of the sliding window. default=10000

=item B<--step-size>

the size of the sliding window steps. default=1000


=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

 a sorted pileup file; Sorting has to be done by (i) ascii reference-id and (ii) numerical position in the ref-seq
 This can be achieved using the Unix-command:
 
 sort -k 1,1 -k 2,2n unsorted.pileup > sorted.pileup
 
 example of correct sorting:
 2L 2LHet 2R 2RHet 3L 3LHet 3R 3RHet 4 U_minus_mitoch Uextra X XHet YHet dmel_mitochondrion_genome gi|42519920|ref|NC_002978.6
 
=head2 Output

 this program generate 2 types of output
 1) chromosome wise output
 
 2L	84.00	37.55	62.45
 X	94.90	39.82	60.18
 
 col1: reference contig (chromosome)
 col2: Average coverage for contig (chromosome)
 col3: GC content in contig (chromosome)
 col4: AT content in contig (chromosome)

 2) sliding window wise output
 
 2L	5000	0.51	na	na	na
 2L	10000	1.00	85.32	37.48	62.52
 2L	15000	0.99	78.98	35.68	64.32
 
 col1: reference contig (chromosome)
 col2: mean position of the sliding window
 col3: fraction of the window which is covered by at least one read; if a fraction of a window lower than defined by the C<--min-covered-fraction> parameter is covered, average coverage, GC and AT content values will not be calculated for that window; the average across a chromosome is independent of the fraction covered
 col4: Average coverage of the sliding window
 col5: GC content of the sliding window
 col6: AT content of the sliding window
 
=head2 Description

 
=head1 AUTHORS

Ram vinay pandey

=cut
