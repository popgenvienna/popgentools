{
    package ParseSam;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA    = qw(Exporter);
    our @EXPORT = qw(get_basic_samparser get_te_samparser get_basic_samformater);
    
    
#@SQ     SN:X    LN:20661737
#HWUSI-EAS300R:2:1:2:807#0       99      2R      9302642 60      73M     =       9302816 220     TAA
#AGTTAAATAATGTACAGTTTAATAACTTAATAACAGTGCAGTTAATTTGAGTTTAAGCGAATCAAACAAG       BBBBBBCBBACBBCBBBBBBAB
#CCBBCBBBCBBBBBBBB>BBBBABBABBBBBBABBB@BBCAABBBABAB>@       XT:A:U  NM:i:0  SM:i:37 AM:i:37 X0:i:1  X
#1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:73
    
    sub get_basic_samformater
    {
        return sub
        {
            my $sam=shift;
            my @entry=($sam->{readid},$sam->{flag},$sam->{chr},$sam->{start},$sam->{mq},$sam->{cigar},$sam->{chrmate},$sam->{posmate},$sam->{distance},$sam->{seq},$sam->{qual},$sam->{appendix});
            return join("\t",@entry); 
        }
    }



    sub get_basic_samparser
    {
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
            
            ## temporary adaption
            #$entry->{chrmate}=~s/_.+//;
            #$entry->{chr}=~s/_.+//;
          return $entry;
         }
    }
    
    sub get_te_samparser
    {
        my $te_resolver=shift;
        my $bp=get_basic_samparser();

        
        # readid, flag, chr, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
        # start, start_s, end, end_s, te_ins_read, te_ins_direction, pp
        return sub
        {
            my $l=shift;
            my $s=$bp->($l);
            
            return $s if $s->{flag} & 0x004 or $s->{flag} & 0x008;
            my $start=$s->{start};
            my $cigar=$s->{cigar};
            
            my ($leng,$startoverhang,$endoverhang)=_getAlignmentLength($cigar);
            my ($end,$startshade, $endshade)=($start+$leng-1, $start-$startoverhang, $start+$leng+$endoverhang-1);
            $s->{end}=$end;
            $s->{start_s}=$startshade;
            $s->{end_s}=$endshade;
            
            my $teinsread=0;
            my $direction="-";
            my $chr=$s->{chr};
            my $chrmate=$s->{chrmate};
            die "both reads need to be mapped"if $chr eq "*" or $chrmate eq "*";
            $s->{pp}=0;
            
            
            if($chrmate eq "=")
            {
                unless($te_resolver->($chr))
                {
                    my $flag=$s->{flag};
                    if((($flag & 0x0010) != ($flag & 0x0020)))
                    {
                        # check if the forward inserted read is before the reverse inserted read
                        my($sfwd,$srev)=($flag&0x0010)?($s->{posmate},$s->{start}):($s->{start},$s->{posmate});
                        $s->{pp}=1 if $sfwd<$srev;
                    }
                }
            }
            else
            {
                my $chrmate_exist=$te_resolver->($chrmate);
                my $chr_exist=$te_resolver->($chr);
                
                if($chrmate_exist and not $chr_exist)
                {
                   # the mate is on a TE; this read is not on a TE (i.e: it is on a proper contig)
                   $teinsread=1;
                   $direction=$s->{flag} & 0x0010?"R":"F";  
                }
            }

            $s->{te_ins_read}=$teinsread;
            $s->{te_ins_direction}=$direction;
            return $s;
        }
        
    }
    
    
    sub _getAlignmentLength
    {
        my $cigar=shift;
        die "no proper cigar provided" unless $cigar;
        die "Can not deal with hard clipping or padding"  if $cigar=~/[HP]/;
        my (@e)=$cigar=~/(\d+[MSDIN])/g;
        
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


    

}


1;