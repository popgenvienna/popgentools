{
    package ParseSam;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA    = qw(Exporter);
    our @EXPORT = qw(parseSam calculate_sam_statistics);
    
    
#@SQ     SN:X    LN:20661737
#HWUSI-EAS300R:2:1:2:807#0       99      2R      9302642 60      73M     =       9302816 220     TAA
#AGTTAAATAATGTACAGTTTAATAACTTAATAACAGTGCAGTTAATTTGAGTTTAAGCGAATCAAACAAG       BBBBBBCBBACBBCBBBBBBAB
#CCBBCBBBCBBBBBBBB>BBBBABBABBBBBBABBB@BBCAABBBABAB>@       XT:A:U  NM:i:0  SM:i:37 AM:i:37 X0:i:1  X
#1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:73




    sub parseSam
    {
        my $line=shift;
        my @a=split /\t/,$line;
        
        my $entry=
        {
            readid=>$a[0],
            flag=>$a[1],
            chr=>$a[2],
            pos=>$a[3],
            mq=>$a[4],
            cigar=>$a[5],
            chrmate=>$a[6],
            posmate=>$a[7],
            distance=>$a[8],
            seq=>$a[9],
            qual=>$a[10],
            p=>$a[1] & 0x0001, # read is paired in sequence
            P=>$a[1] & 0x0002, # mapped in a proper pair
            r=>$a[1] & 0x0010?"R":"F", # strand of the query (1 for reverse, 0 for forward)
            R=>$a[1] & 0x0020?"R":"F", # strand of the mate (1 for reverse, 0 for forward)
        };
        
        return $entry;
    }
    
    
    sub calculate_sam_statistics
    {
        my $data=shift;
        my $countall=0;
        my($mappedInProperPair,$fucked,$fuckedMateUnmapped,$fuckedSamePosition)=(0,0,0,0,0);
        my($fuckedSameStrand,$fuckedMatesOverlap,$fuckedImpropableDistance,$fuckedOtherContig,$countSingleEnd)=(0,0,0,0,0);
        
        return
        {
            mall=>"na",
            mpp=>"na",
            mse=>"na",
            f_mipp=>"na",
            f_mum=>"na",
            f_sp=>"na",
            f_ss=>"na",
            f_mo=>"na",
            f_id=>"na",
            f_oc=>"na"        
        } unless @$data;
        
        foreach my $read (@$data)
        {
            
            my $flag        =$read->{flag};
            my $matechr     =$read->{chrmate};
            my $matedist    =$read->{distance};
            my $matepos     =$read->{posmate};
            my $refpos      =$read->{pos};
            my $readlength  =length($read->{seq});
            
            # discard crap
            next if $flag & 0x0004;
            $countall++;
            
            if($flag & 0x0001)
            {
                # paired in sequence?
                if($flag & 0x0002)
                {
                    # proper paired in mapping
                    $mappedInProperPair++;
                }
                else 
                {
                    # not-proper pair in mapping
                    $fucked++;
                    
                    if($flag & 0x0008) # only one read is mapped
                    {
                        # mate unmapped
                        $fuckedMateUnmapped++;
                    }
                    else 
                    {
                        # a.) not-proper pair in mapping
                        # b.) mate is mapped
                        if($matechr eq "=") # same contig
                        {
                            # mate is mapped to the same chromosome
                            if($matepos==$refpos)
                            {
                                # the mate has the same position -> identical
                                $fuckedSamePosition++;
                            }
                            else
                            {
                                # mate does not have the same position   
                                my $stranda=$flag & 0x0010?"R":"F";
                                my $strandb=$flag & 0x0020?"R":"F";
    
                                if($stranda eq $strandb)
                                {
                                    # mate is on the same strand
                                    $fuckedSameStrand++;
                                }
                                else
                                {
                                    # mate is on a different strand
                                    my $absdist=abs($matedist);
                                    if ($absdist<$readlength)
                                    {
                                        # mate is overlapping
                                        $fuckedMatesOverlap++;
                                    }
                                    else
                                    {
                                        # mate has a weired distance
                                        $fuckedImpropableDistance++;
                                    }
                                }
                                
                            }
                            
                        }
                        else  
                        {
                            # the mate is mapped to a different contig
                            $fuckedOtherContig++;
                        }
                    }
    
                }
    
            }
            else
            {
                $countSingleEnd++;
            }
        }
        
        #my($mappedInProperPair,$fucked,$fuckedMateUnmapped,$fuckedSamePosition)=(0,0,0,0,0);
        #my($fuckedSameStrand,$fuckedMatesOverlap,$fuckedImpropableDistance,$fuckedOtherContig,$countSingleEnd)=(0,0,0,0,0);
        return
        {
            mall=>$countall,   # mapped all
            mpp=>$mappedInProperPair,
            mse=>$countSingleEnd, # single ended reads mapped as single end
            f_mipp=>$fucked, # fucked: mapped in improper pair
            f_mum=>$fuckedMateUnmapped,
            f_sp=>$fuckedSamePosition,
            f_ss=>$fuckedSameStrand,
            f_mo=>$fuckedMatesOverlap,
            f_id=>$fuckedImpropableDistance,
            f_oc=>$fuckedOtherContig        
        };
    }

}


1;