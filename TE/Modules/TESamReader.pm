{
    package TESamReader;
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";

    
    
    sub new
    {
        my $class=shift;
        my $file=shift;
        my $sp=shift;
        my $makenoise=shift;
        
        open my $fh,"<",$file or die "Could not open file handle";
        
        
        my $tobless={
            noise   =>$makenoise,
            counter =>0,
            file    =>$file,
            pair    =>{},
            fh      =>$fh,
            samparser=>$sp,
            buffer  =>[],
        };
        my $self=bless $tobless, __PACKAGE__;
        
        # spool forward; ie. skip the header
        while(1)
        {
            my $l=$self->_nextline();
            unless($l=~/^@/)
            {
                $self->_bufferline($l);
                last;
            }
        }
        return $self;
        
    }
    
    
    sub next
    {
        my $self=shift;

        my $sp = $self->{samparser};
        my $ph=$self->{pair};
        
        while(1)
        {
            my $line=$self->_nextline();
            $self->{counter}++;
            print "Processed $self->{counter} sam entries\n" unless $self->{counter} % $self->{noise};
            last unless $line;
            chomp $line;
            my $s=$sp->($line);
            # readid, flag, chr, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
            # start, start_s, end, end_s, te_ins_read, te_ins_direction, pp
            
            # discard unmapped or pairunmapped
            next if $s->{flag} &0x004;
            next if $s->{flag} &0x008;
            
            if($s->{te_ins_read})
            {
                return
                {
                  mode=>"pre",
                  r1=>$s
                };
            }
            elsif($s->{pp})
            {
                my $rid=$s->{readid};
                if(exists($ph->{$rid}))
                {
                    
                    my ($r1,$r2)=($ph->{$rid},$s);
                    
                    # r1 should always have a lower start position as r2;
                    # thus r1 must always be the fwd insertion and r2 always the reverse insertion read
                    ($r1,$r2)=($r2,$r1) if($r1->{start_s} > $r2->{start_s});
                    delete($ph->{$rid}); # delete the hash entry -> avoid memory overkill
                    return
                    {
                        mode=>"abs",
                        r1=>$r1,
                        r2=>$r2
                    };
                }
                else
                {
                    $ph->{$rid}=$s;
                }
            }
        }
        return undef;
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
    
    sub close
    {
        my $self=shift;
        my $fh=$self->{fh};
        close $fh;
    }
}

1;