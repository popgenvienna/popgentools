{
    package Test_SamSlider;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use SamSlider;
    use Test;
    
    require Exporter;
    our @ISA    = qw(Exporter);
    our @EXPORT = qw(test_SamSlider);
    
#HWUSI-EAS300R:2:1:2:807#0       99      2R      9302642 60      73M     =       9302816 220     TAA
#AGTTAAATAATGTACAGTTTAATAACTTAATAACAGTGCAGTTAATTTGAGTTTAAGCGAATCAAACAAG       BBBBBBCBBACBBCBBBBBBAB
#CCBBCBBBCBBBBBBBB>BBBBABBABBBBBBABBB@BBCAABBBABAB>@       XT:A:U  NM:i:0  SM:i:37 AM:i:37 X0:i:1  X
#1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:73
    sub test_SamSlider
    {
        my $ss; # SamSlider;
        my $sas; # SamString
        my $w;      # window
        $sas=
        "R1\t16\t2R\t5\t30\t73M\t=\t10\t5\tAAAA\tbbbb\tsheit1\n".
        "R2\t0\t2R\t6\t31\t75M\t=\t11\t20\tAACC\tbbaa\tsheit2\n".
        "R3\t0\t2R\t7\t34\t73M\t=\t10\t13\tGACC\tabaa\tsheit3\n".
        "R4\t0\t2R\t8\t29\t73M\t=\t10\t5\tGGGG\tccaa\tsheit4\n".
        "R5\t16\t3R\t3\t40\t73M\t=\t10\t12\tAACC\tbbaa\tsheit5\n".
        "R6\t0\t4R\t4\t31\t73M\t=\t10\t5\tAACC\tbbaa\tsheit2\n".
        "R7\t16\t4R\t5\t31\t73M\t=\t10\t5\tAACC\tbbaa\tsheit2\n".
        "R8\t0\t4R\t6\t31\t73M\t=\t10\t5\tAACC\tbbaa\tsheit2\n".
        "R9\t0\t4R\t7\t31\t73M\t=\t10\t5\tAACC\tbbaa\tsheit2\n";
        $ss=_getSamSliderFromString($sas,3,1);
        $w=$ss->nextWindow();
        is($w->{chr},"2R","test SamSlider; reference ok");
        is($w->{start},"0","test SamSlider; start position ok");
        is($w->{end},"3","test SamSlider; end position ok");
        is($w->{window},"3","test SamSlider; window size ok");
        is($w->{count},"0","test SamSlider; data count ok");
        
        $w=$ss->nextWindow();
        $w=$ss->nextWindow();
        $w=$ss->nextWindow();
        is($w->{chr},"2R","test SamSlider; reference ok");
        is($w->{start},"3","test SamSlider; start position ok");
        is($w->{end},"6","test SamSlider; end position ok");
        is($w->{window},"3","test SamSlider; window size ok");
        is($w->{count},"2","test SamSlider; data count ok");
        is($w->{data}[0]{pos},"5","test SamSlider; data position is ok");
        is($w->{data}[1]{pos},"6","test SamSlider; data position is ok");
        is($w->{data}[0]{chr},"2R","test SamSlider; chromosome is ok");
        is($w->{data}[0]{cigar},"73M","test SamSlider; cigar is ok");
        is($w->{data}[0]{distance},"5","test SamSlider; distance position is ok");
        is($w->{data}[0]{flag},"16","test SamSlider; binary flag is ok");
        is($w->{data}[0]{mq},"30","test SamSlider; mapping quality is ok");
        is($w->{data}[0]{readid},"R1","test SamSlider; readid is ok");
        is($w->{data}[0]{r},"R","test SamSlider; strand of the query is ok");
        is($w->{data}[0]{posmate},"10","test SamSlider; position of the mate is ok");
        is($w->{data}[0]{seq},"AAAA","test SamSlider; sequence of the read is ok");
        is($w->{data}[0]{qual},"bbbb","test SamSlider; quality sequence is ok");
        
        $w=$ss->nextWindow();
        $w=$ss->nextWindow();
        is($w->{chr},"2R","test SamSlider; reference ok");
        is($w->{start},"5","test SamSlider; start position ok");
        is($w->{end},"8","test SamSlider; end position ok");
        is($w->{window},"3","test SamSlider; window size ok");
        is($w->{count},"3","test SamSlider; data count ok");
        is($w->{data}[0]{pos},"6","test SamSlider; data position is ok");
        is($w->{data}[1]{pos},"7","test SamSlider; data position is ok");
        is($w->{data}[2]{pos},"8","test SamSlider; data position is ok");
        is($w->{data}[0]{chr},"2R","test SamSlider; chromosome is ok");
        is($w->{data}[0]{cigar},"75M","test SamSlider; cigar is ok");
        is($w->{data}[0]{distance},"20","test SamSlider; distance position is ok");
        is($w->{data}[0]{flag},"0","test SamSlider; binary flag is ok");
        is($w->{data}[0]{mq},"31","test SamSlider; mapping quality is ok");
        is($w->{data}[0]{readid},"R2","test SamSlider; readid is ok");
        is($w->{data}[0]{r},"F","test SamSlider; strand of the query is ok");
        is($w->{data}[0]{posmate},"11","test SamSlider; position of the mate is ok");
        is($w->{data}[0]{seq},"AACC","test SamSlider; sequence of the read is ok");
        is($w->{data}[0]{qual},"bbaa","test SamSlider; quality sequence is ok");
       
        # R4\t0\t2R\t8\t29\t73M\t=\t10\t5\tGGGG\tccaa\tsheit4\n".
        is($w->{data}[2]{chr},"2R","test SamSlider; chromosome is ok");
        is($w->{data}[2]{cigar},"73M","test SamSlider; cigar is ok");
        is($w->{data}[2]{distance},"5","test SamSlider; distance position is ok");
        is($w->{data}[2]{flag},"0","test SamSlider; binary flag is ok");
        is($w->{data}[2]{mq},"29","test SamSlider; mapping quality is ok");
        is($w->{data}[2]{readid},"R4","test SamSlider; readid is ok");
        is($w->{data}[2]{r},"F","test SamSlider; strand of the query is ok");
        is($w->{data}[2]{posmate},"10","test SamSlider; position of the mate is ok");
        is($w->{data}[2]{seq},"GGGG","test SamSlider; sequence of the read is ok");
        is($w->{data}[2]{qual},"ccaa","test SamSlider; quality sequence is ok");
        
        # "R5\t16\t3R\t3\t40\t73M\t=\t10\t12\tAACC\tbbaa\tsheit5\n"
        $w=$ss->nextWindow();
        is($w->{chr},"3R","test SamSlider; reference ok");
        is($w->{start},"0","test SamSlider; start position ok");
        is($w->{end},"3","test SamSlider; end position ok");
        is($w->{window},"3","test SamSlider; window size ok");
        is($w->{count},"1","test SamSlider; data count ok");
        is($w->{data}[0]{chr},"3R","test SamSlider; chromosome is ok");
        is($w->{data}[0]{cigar},"73M","test SamSlider; cigar is ok");
        is($w->{data}[0]{distance},"12","test SamSlider; distance position is ok");
        is($w->{data}[0]{flag},"16","test SamSlider; binary flag is ok");
        is($w->{data}[0]{mq},"40","test SamSlider; mapping quality is ok");
        is($w->{data}[0]{readid},"R5","test SamSlider; readid is ok");
        is($w->{data}[0]{r},"R","test SamSlider; strand of the query is ok");
        is($w->{data}[0]{posmate},"10","test SamSlider; position of the mate is ok");
        is($w->{data}[0]{seq},"AACC","test SamSlider; sequence of the read is ok");
        is($w->{data}[0]{qual},"bbaa","test SamSlider; quality sequence is ok");
        
        $w=$ss->nextWindow();
        is($w->{chr},"4R","test SamSlider; reference ok");
        is($w->{start},"0","test SamSlider; start position ok");
        is($w->{end},"3","test SamSlider; end position ok");
        is($w->{window},"3","test SamSlider; window size ok");
        is($w->{count},"0","test SamSlider; data count ok");
        $w=$ss->nextWindow();
        is($w->{chr},"4R","test SamSlider; reference ok");
        is($w->{start},"1","test SamSlider; start position ok");
        is($w->{end},"4","test SamSlider; end position ok");
        is($w->{window},"3","test SamSlider; window size ok");
        is($w->{count},"1","test SamSlider; data count ok");
        $w=$ss->nextWindow();
        $w=$ss->nextWindow();
        $w=$ss->nextWindow();
        is($w->{chr},"4R","test SamSlider; reference ok");
        is($w->{start},"4","test SamSlider; start position ok");
        is($w->{end},"7","test SamSlider; end position ok");
        is($w->{window},"3","test SamSlider; window size ok");
        is($w->{count},"3","test SamSlider; data count ok");
        
        $w=$ss->nextWindow();
        not_exists($w,"test SamSlider; end of file identified correctly");
        my $bla=0;
    }
    
    sub _getSamSliderFromString
    {
        my $str=shift;
        my $window=shift;
        my $step=shift;

        open my $ofh,"<",\$str or die "could not open string filehandle";
        my $cr=bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            fh=>$ofh,
            curwin=>[],
            buffer=>[]
        },"SamSlider";
        return $cr;
    }
    
    
    
}