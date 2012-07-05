{
    package SamSlider;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use ParseSam;
    
    
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

# important for sam field
# chr pos 

    
    sub nextWindow
    {
        my $self=shift;
        
        #get the current window, and the current chromosome
        my $curwin=$self->{curwin};
        
        my $curChr="";
        $curChr=$curwin->[0]{chr} if @$curwin;
        
        my $resetchr=0;
        
        # empty unnecessary entries
        EMPTY: while(@$curwin)
        {
            my $e=shift @$curwin;
            if($e->{pos} > $self->{lower})
            {
                unshift @$curwin, $e;
                last EMPTY;
            }
            
        }
        
        # fill with novel entries
        my $line;
        FILL:while($line=$self->_nextline)
        {
            # skip the SAM header
            next FILL if $line=~/^@/;
            
            my $e=parseSam($line);
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
        
        
        my $toret=$self->_annotateWindow($curwin,$curChr,$self->{lower},$self->{upper},$self->{window});
        
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
    
    sub _annotateWindow
    {
        my $self=shift;
        my $curwin=shift;
        my $chr=shift;
        my $start=shift;
        my $end=shift;
        my $windowsize=shift;
        
        my $win;
        my $count=0;
        foreach my $e(@$curwin)
        {
            $count++;
            die "chromosomes are not identical" if $e->{chr} ne $chr;
        }
        
        $win={
          chr=>$chr,
          start=>$start,
          end=>$end,
          count=>$count,
          data=>$curwin,
          middle=>(($end+1+$start)/2),
          window=>$windowsize
        };
        return $win;
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


1;