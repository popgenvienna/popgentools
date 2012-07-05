{
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use List::Util qw/max min sum/;


    my $help=0;
    my $test=0;
    my $threshold=0.0;
    my $minScore=0.0;
    my $input="";
    my $output="";
    my $verbose=1;
    my $pd_low=0;

    
    GetOptions(
        "input=s"               =>\$input,
        "output=s"              =>\$output,
        "threshold=f"           =>\$threshold,
        "min-score=f"           =>\$minScore,
        "quit"                  =>sub{$verbose=0},
        "test"                  =>\$test,
        "help"                  =>\$help
    ) or pod2usage(-msg=>"Options not valid $!",-verbose=>1);
    
    pod2usage(-verbose=>2) if $help;
    TestTrim::runTrimTests() if $test;
    pod2usage(-msg=>"A input file has to be provided", -verbose=>1) unless -e $input;
    pod2usage(-msg=>"A output file has to be provided", -verbose=>1) unless  $output;



    
    open my $ifh,"<",$input or die "Could not open input file";
    open my $ofh,">",$output or die "Could not open output file";
    
    my $scorecalc=get_score_calculator($pd_low,$threshold);
    my $hsslider=HighScoreSlider->new($input,$scorecalc);
    
    
    while(my $hsp=$hsslider->next())
    {
        next unless ($hsp->{hs});
        next if $hsp->{hs} < $minScore;
        my @vals=sort {$b->{val}<=>$a->{val}} @{$hsp->{val}};
        my $entries =@vals;
        my $max=$vals[0];
        
        my $extreme=$pd_low ? min(@vals) : max(@vals);
        my $av=sum(map {$_->{val}} @vals)/$entries;
        
        my $hs=sprintf "%.6f",$hsp->{hs};
        $av=sprintf "%.6f",$av;
        $extreme=sprintf "%.6f",$extreme;

        
        
        
        print $ofh "$hsp->{chr}\t$hsp->{start}\t$hsp->{end}\t$entries\t$hs\t$av\t$max->{val}\t$max->{pos}\n";
    }
    
exit;    
}






    

sub get_score_calculator
{
    my $pdlow=shift;
    my $threshold=shift;
    
    
    # Low
    # t=0.1 v=0.05 -> s=0.05
    # t=0.1 v=0.2  -> s=-0.1
    # t=-1  v=-2.2  -> s=1.2
    # t=-1  v=0.3   -> s=-1.3
    # score =-(v-t)
    
    # High
    # t=0.5 v=0.7  -> s=0.2
    # t=0.5 v=0.2  -> s=-0.3
    # t=-1  v=0.5  -> s=1.5
    # t=-1  v=-2.2 -> s=-1.2
    # score = v-t
    
    
    return sub
    {
        my $val=shift;
        my $score=$val-$threshold;
        $score= (-1)* $score if $pdlow;
        return $score;
    }
}

{
    package HighScoreSlider;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    
    
    sub new
    {
        #($input,$stepsize,$scorecalc)
        my $class       =shift;
        my $file        =shift;
        my $scorecalc   =shift;

        open my $fh,"<$file" or die "Could not open file handle";        
        return bless {
            file        => $file,
            fh          => $fh,
            scorecalc   => $scorecalc,
            buffer=>[]
        },__PACKAGE__;
    }
    

    sub next
    {
        my $self=shift;
        while(my $hsc=$self->_nextWin())
        {
            return $hsc if $hsc->{hs}>0;
        }
        return undef;
    }
    
    sub _nextWin
    {
        my $self=shift;

    
        # fill with novel entries
        my $line;
        
        my $curChr;
        my $highscorestart=-1;
        my $highscoreend=-1;
        my $highscore=0;
        my $runningscore=0;
        my $lastpos="na";
        my $highscorevals=[];
        my $counter=0;
        my $highscorecounter=0;
        
        while($line=$self->_nextline)
        {
            chomp $line;
            #parse line
            my $e=$self->_parseVariance($line);
            my $chr=$e->{chr};
            my $pos=$e->{pos};
            my $val=$e->{val};
            next if $val eq "na";
            # set the controls
            $curChr = $chr unless $curChr;

            # check the abortion conditions;
            # 1..chromosome change
            # 2..position is not increment of the step size
            # 3..value equals to "na"
            if($curChr ne $e->{chr})
            {
                ######################################################################################################
                $self->_bufferline($line);
                return $self->_annotateHighscore($curChr, $highscorestart, $highscoreend, $highscore, $highscorevals, $highscorecounter);
            }
            
            # score lower than zero -> abort
            my $s=$self->{scorecalc}->($val);
            $runningscore+=$s;
            ###########################################################################################################################
            return $self->_annotateHighscore($curChr, $highscorestart, $highscoreend, $highscore, $highscorevals, $highscorecounter) if ($runningscore<0);
            
            
            
            
            # everything passing here may be potentiall be contributing to the new highscore
            push @$highscorevals,$e; #only the first fraction actually contributing to the highscore must be considered
            
            $counter++;
            # score higher than zero
            if($runningscore>$highscore and $highscorestart==-1)
            {
                $highscorevals=[$e];
                $highscorestart=$pos;
                $highscore=$runningscore;
                $highscoreend=$pos;
                $counter=1;
                $highscorecounter=1;
            }
            elsif($runningscore>$highscore)
            {
                $highscoreend=$pos;
                $highscore=$runningscore;
                $highscorecounter=$counter;
            }
        }
        
        
        if($highscore)
        {
            #####################################################################################################
           return $self->_annotateHighscore($curChr, $highscorestart, $highscoreend, $highscore, $highscorevals,$highscorecounter);
        }
        else
        {
            return undef;
        }
        
        
    }
    
    
    sub _parseVariance
    {
        my $self=shift;
        my $line=shift;
        #2L      27      28      fst_single      0.00518329      0.00257084      0.00281954      0.183731
        #2L      141     142     fst_single      0.00788214      0.00775194      0.00222717      0.544504
        #2L      334     335     fst_single      0.00699301      0.00028161      0.00671141      -0.110606
        #2L      493     494     fst_single      0.01233904      0.00239368      0.00488322      0.327132
    
        my ($chr,$pos,undef,undef,undef,undef,undef,$val)=split /\t/,$line;
        die "Position is not numeric $pos" unless $pos=~/^\d+$/;
        return {
            chr=>$chr,
            pos=>$pos,
            val=>$val
        };
        
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
    
    
    sub _annotateHighscore
    {
        my $self=shift;
        my $chr=shift;
        # _annotateHighscore($highscorestart, $highscoreend, $highscore, $highscorevals);
        my $start=shift;
        my $end=shift;
        my $highscore=shift;
        my $vals=shift;
        my $entries=shift;
        unless($highscore>0)
        {
            return
            {
                hs=>0,
                start=>0,
                end=>0
            };
        }
        
        die "impossible number of entries" unless $entries;
        
        
        die "must not happen" if $end<$start;
        

        my $realstart=$start;
        my $realend=$end;
        
        #extract subhighscores
        my @vvals=@$vals;
        my @hsvals=@vvals[0..$entries-1];
        die "not allowed" unless @hsvals ==$entries;
        
        
        return
        {
            chr=>$chr,
            start=>$realstart,
            end=>$realend,
            hs=>$highscore,
            val=>\@hsvals,
        };
    }
    
    
    

}




=head1 NAME

Identify-variance-peaks.pl - Identify peaks in the output of Variance-sliding using a Mott algorithm

=head1 SYNOPSIS

 perl Identify-variance-peaks --input input.pi --output peaks --step-size 1000 --threshold 0.05 --min-score 0.4 
 
=head1 OPTIONS

=over 4

=item B<--input>

The input file; A file produced by Variance-slider.pl. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--step-size>

The step size which has been used to run Variance-slider; Mandatory parameter

=item B<--threshold>

threshold for calculating the score; score = variance-threshold; default=0

=item B<--min-score>

Minimum score for a variance peaks; Variance peaks with a score smaller than this will be discared; default=0

=item B<--high-peak>

flag; the script per default identifies peaks of low variance (e.g.: pi, Tajimas'D). For example: the more negative Tajima's D the more pronounced is the signature of positive slection;
However, if someone is interested instead in balancing selection he would like to see peaks of high Tajima's D; this flag allows to swith from the default (low) to high peak scoreing. default=off


=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS

=head2 INPUT

The output of Variance Sliding

 2L      25000   7       0.829   0.000770972
 2L      35000   16      0.978   0.000779493
 2L      45000   34      0.757   0.002023476
 2L      55000   19      0.733   0.001450490
 2L      65000   28      0.960   0.001053402



=head2 OUTPUT

 2L      1390001 1420000 0.003638        0.007790        0.006362
 2L      1460001 1490000 0.004169        0.007957        0.005831
 2L      1500001 1540000 0.005756        0.006734        0.004244
 2L      1560001 1570000 0.003235        0.006765        0.006765
 
 col 1: contig
 col 2: start
 col 3: end
 col 4: score for the peak
 col 5: average measure for the peak (pi, theta, Tajima's D)
 col 6: the most pronounced measure for the peak (pi, theta, Tajima's D)

=cut