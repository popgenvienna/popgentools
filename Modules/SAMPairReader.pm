{
    package SAMPairReader;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";

    
    sub _get_basic_samparser
    {
        my $class=shift;
        return sub
        {
            my $line=shift;
            my @a=split /\t/,$line;
            
            # readid, flag, chr, start, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
            my $entry=
            {
                readid  =>$a[0],
                flag    =>$a[1],
                chr     =>$a[2],
                start   =>$a[3],
                mq      =>$a[4],
                cigar   =>$a[5],
                chrmate =>$a[6],
                posmate =>$a[7],
                distance=>$a[8],
                seq     =>$a[9],
                qual    =>$a[10],
                appendix=>$a[11]
          };
            
            my ($leng,$startoverhang,$endoverhang)=(1,0,0);
            ($leng,$startoverhang,$endoverhang)=_getAlignmentLength($entry->{cigar}) if $entry->{cigar} ne "*";
            my ($end,$startshade, $endshade)=($entry->{start}+$leng-1, $entry->{start}-$startoverhang, $entry->{start}+$leng+$endoverhang-1);
            $entry->{end}=$end;
            $entry->{start_s}=$startshade;
            $entry->{end_s}=$endshade;
          return $entry;
         }
    }
    
    sub _getAlignmentLength
    {
        my $cigar=shift;
        die "no proper cigar provided" unless $cigar;
        die "Can not deal with hard clipping or padding"  if $cigar=~/[HP]/;
        my (@e)=$cigar=~/(\d+[MSDIN])/g;
        die "$cigar not appropriate" unless $cigar=~m/^(\d+[MSDIN])+$/;
        
        # get rid of the first shaded thing
        #HWUSI-EAS300R:7:1:8:239#0       147     3R      10371697        24      3S70M1S =       10371524        -173    TCGGCCATAAGTTCACTTGAGAGCAGTGACATAACCAGGGTTTTCCAGGGCGACAGGTTAGCTGAAGGTAGGTN
        #      BBBBBBBBBYUWSNKNYNKPMWWWYWWWNWYTSQWTNWYYTTUWWWWWWWWLPNUQUUYNWTTWWWYYWWWWOD      AS:i:58
        
        my $startoverhang=0;
        if($e[0]=~/^(\d+)S$/)
        {
            $startoverhang=$1;
            shift @e;
        }
        
        my $endoverhang=0;
        if($e[-1]=~/^(\d+)S$/)
        {
            $endoverhang=$1;
            pop @e;
        }
        
        my $alrefleng=0; # alignment length with respect to the reference
        foreach my $p (@e)
        {
            if($p=~/^(\d+)[M]$/)
            {
                $alrefleng+=$1;
            }
            elsif($p=~/^(\d+)I/)
            {
                # insertions; insertions in the reads do not have a consequence on the position in the reference -> ignore
            }
            elsif($p=~/^(\d+)[ND]/)
            {
                $alrefleng+=$1;
            }
            else
            {
                die "cigar entry $p not supported";
            }
        }
        return ($alrefleng,$startoverhang,$endoverhang);
    }
    
    sub new
    {
        my $class=shift;
        my $file=shift;
        my $makenoise=shift;
        
        open my $fh,"<",$file or die "Could not open file handle";
        
        
        my $tobless={
            noise   =>$makenoise,
            counter =>0,
            file    =>$file,
            pair    =>{},
            fh      =>$fh,
            buffer  =>[],
        };
        my $self=bless $tobless, __PACKAGE__;
        $self->{samparser}=$self->_get_basic_samparser();
        
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
        my $ph = $self->{pair};
        
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
            

            my $rid=$s->{readid};
            if(exists($ph->{$rid}))
            {
                    
                 my ($r1,$r2)=($ph->{$rid},$s);
                    
                # the first read should always map to the forward strand     
                ($r1,$r2)=($r2,$r1) if($r1->{flag} & 0x0010);
                delete($ph->{$rid}); # delete the hash entry -> avoid memory overkill
                return [$r1,$r2];

            }
            else
            {
                $ph->{$rid}=$s;
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
