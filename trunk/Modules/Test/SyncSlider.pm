{
    package Test::SyncSlider;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use SyncSlider;
    use Synchronized;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT=qw(run_SyncSliderTests);
    our @EXPORT_OK = qw();
    
    
    sub run_SyncSliderTests
    {
        testSyncSlider();
    }
    
    
    sub _getSyncSliderForString #str, window, step, mincount, mincov, maxcov
    {
        my $str=shift;
        my $window=shift;
        my $step=shift;
        my $mincount=shift;
        my $mincov=shift;
        my $maxcov=shift;
        $maxcov ||=1000000;
        my $sp=get_sumsnp_synparser($mincount,$mincov,$maxcov);

        open my $ofh,"<",\$str or die "could not open string filehandle";
        my $cr=bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            fh=>$ofh,
            sp=>$sp,
            curwin=>[],
            buffer=>[]
        },"SyncSlider";
        return $cr;
    }
    
    sub testSyncSlider
    {
        my $teststr=
        "chr1\t1\tN\t0:0:6:0:0:0\t0:0:0:7:0:0\t0:0:8:0:0:0\n".
        "chr1\t2\tN\t0:4:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n".
        "chr1\t3\tN\t0:0:3:0:0:0\t0:0:0:8:0:0\t0:0:13:0:0:0\n".
        "chr1\t4\tN\t0:0:5:0:0:0\t0:0:4:0:0:0\t0:0:6:0:0:0\n".
        "chr1\t5\tN\t0:0:8:0:0:0\t0:0:0:8:0:0\t0:0:8:0:0:0\n".
        "chr2\t3\tN\t0:0:11:0:0:0\t0:0:0:8:0:0\t0:0:20:0:0:0\n".
        "chr3\t1\tN\t0:0:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n".
        "chr3\t2\tN\t0:0:1:0:0:0\t0:0:0:5:0:0\t0:0:8:0:0:0\n".
        "chr4\t4\tN\t0:0:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n";
        my $bpsl=_getSyncSliderForString($teststr,3,1,2,4);
        
        my $w=$bpsl->nextWindow();
        
        is($w->{chr},"chr1","test Synclider, correct chromosome");
        is($w->{count_covered},2,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},2,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        is($w->{data}[0]{ispuresnp},1,"test SyncSlider, correct identification of a pure SNP at the first position");
        is($w->{data}[1]{ispuresnp},1,"test SyncSlider, correct identification of a pure SNP at the second position");
        is($w->{data}[2]{ispuresnp},0,"test SyncSlider, correct, there is no pure SNP at position two");
        is($w->{data}[0]{issnp},1,"test SyncSlider, correct identification of a SNP at the first position");
        is($w->{data}[1]{issnp},1,"test SyncSlider, correct identification of a SNP at the second position");
        is($w->{data}[2]{issnp},0,"test SyncSlider, correct, there is no SNP at position two");
        is($w->{data}[0]{iscov},1,"test SyncSlider, first position is sufficiently covered");
        is($w->{data}[1]{iscov},1,"test SyncSlider, second position is sufficiently covered for a SNP");
        is($w->{data}[2]{iscov},0,"test SyncSlider, third position is not sufficiently covered for a SNP -> no SNP possible at this position");
        is($w->{data}[0]{pos},1,"test SyncSlider, position is correct");
        is($w->{data}[1]{pos},2,"test SyncSlider, position is correct");
        is($w->{data}[2]{pos},3,"test SyncSlider, position is correct");
        
        is($w->{data}[0]{samples}[0]{eucov},6,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[1]{eucov},7,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[2]{eucov},8,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[0]{totcov},6,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[1]{totcov},7,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[2]{totcov},8,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[0]{index},0,"test SyncSlider, index is correct");
        is($w->{data}[0]{samples}[1]{index},1,"test SyncSlider, index is correct");
        is($w->{data}[0]{samples}[2]{index},2,"test SyncSlider, index is correct");

        is($w->{data}[0]{samples}[0]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[0]{samples}[0]{T},0,"test SyncSlider, count of T is correct");
        is($w->{data}[0]{samples}[0]{C},6,"test SyncSlider, count of C is correct");
        is($w->{data}[0]{samples}[0]{G},0,"test SyncSlider, count of G is correct");
        
        is($w->{data}[1]{samples}[0]{eucov},12,"test SyncSlider, coverage is ok");
        is($w->{data}[1]{samples}[0]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[1]{samples}[0]{T},4,"test SyncSlider, count of T is correct");
        is($w->{data}[1]{samples}[0]{C},8,"test SyncSlider, count of C is correct");
        is($w->{data}[1]{samples}[0]{G},0,"test SyncSlider, count of G is correct");
        
        is($w->{data}[2]{samples}[2]{eucov},13,"test SyncSlider, coverage is ok");
        is($w->{data}[2]{samples}[2]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[2]{samples}[2]{T},0,"test SyncSlider, count of T is correct");
        is($w->{data}[2]{samples}[2]{C},13,"test SyncSlider, count of C is correct");
        is($w->{data}[2]{samples}[2]{G},0,"test SyncSlider, count of G is correct");
        
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr1","test SyncSlider, correct chromosome");
        is($w->{count_covered},2,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test SyncSlider, correct number of snps in region");
        is($w->{start},1,"test SyncSlider, correct start position");
        is($w->{end},4,"test SyncSlider, correct end position");
        is($w->{middle},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        is($w->{data}[0]{pos},2,"test SyncSlider, position is correct");
        is($w->{data}[1]{pos},3,"test SyncSlider, position is correct");
        is($w->{data}[2]{pos},4,"test SyncSlider, position is correct");
        is($w->{data}[0]{ispuresnp},1,"test SyncSlider, correct identification of a SNP at the first position");
        is($w->{data}[1]{ispuresnp},0,"test SyncSlider, correct there is no SNP at the second position");
        is($w->{data}[2]{ispuresnp},0,"test SyncSlider, correct, there is no SNP at position three");
        is($w->{data}[0]{iscov},1,"test SyncSlider, first position is sufficiently covered");
        is($w->{data}[1]{iscov},0,"test SyncSlider, second position is not sufficiently covered for a SNP");
        is($w->{data}[2]{iscov},1,"test SyncSlider, third position is sufficiently covered for a SNP");
        
        is($w->{data}[0]{samples}[0]{eucov},12,"test SyncSlider, coverage is ok");
        is($w->{data}[0]{samples}[0]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[0]{samples}[0]{T},4,"test SyncSlider, count of T is correct");
        is($w->{data}[0]{samples}[0]{C},8,"test SyncSlider, count of C is correct");
        is($w->{data}[0]{samples}[0]{G},0,"test SyncSlider, count of G is correct");
        
        is($w->{data}[1]{samples}[2]{eucov},13,"test SyncSlider, coverage is ok");
        is($w->{data}[1]{samples}[2]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[1]{samples}[2]{T},0,"test SyncSlider, count of T is correct");
        is($w->{data}[1]{samples}[2]{C},13,"test SyncSlider, count of C is correct");
        is($w->{data}[1]{samples}[2]{G},0,"test SyncSlider, count of G is correct");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr1","test SyncSlider, correct chromosome");
        is($w->{count_covered},2,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test SyncSlider, correct number of snps in region");
        is($w->{start},2,"test SyncSlider, correct start position");
        is($w->{end},5,"test SyncSlider, correct end position");
        is($w->{middle},4,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        is($w->{data}[0]{pos},3,"test SyncSlider, position is correct");
        is($w->{data}[1]{pos},4,"test SyncSlider, position is correct");
        is($w->{data}[2]{pos},5,"test SyncSlider, position is correct");
        
        # switch to new chromosome
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr2","test SyncSlider, correct chromosome");
        is($w->{count_covered},1,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{middle},2,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        
        is(scalar(@{$w->{data}}),1,"Correct number of data entries");
        is($w->{data}[0]{pos},3,"test SyncSlider, position is correct");
        is($w->{data}[0]{iscov},1,"test SyncSlider, correct position is sufficiently covered");
        is($w->{data}[0]{ispuresnp},1,"test SyncSlider, correct position is a snp");
        is($w->{data}[0]{samples}[0]{eucov},11,"test SyncSlider, coverage is ok");
        is($w->{data}[0]{samples}[0]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[0]{samples}[0]{T},0,"test SyncSlider, count of T is correct");
        is($w->{data}[0]{samples}[0]{C},11,"test SyncSlider, count of C is correct");
        is($w->{data}[0]{samples}[0]{G},0,"test SyncSlider, count of G is correct");
        is($w->{data}[0]{samples}[2]{eucov},20,"test SyncSlider, coverage is ok");
        is($w->{data}[0]{samples}[2]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[0]{samples}[2]{T},0,"test SyncSlider, count of T is correct");
        is($w->{data}[0]{samples}[2]{C},20,"test SyncSlider, count of C is correct");
        is($w->{data}[0]{samples}[2]{G},0,"test SyncSlider, count of G is correct");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr3","test SyncSlider, correct chromosome");
        is($w->{count_covered},1,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{middle},2,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),2,"Correct number of data entries");
        is($w->{data}[0]{pos},1,"test SyncSlider, position is correct");
        is($w->{data}[1]{pos},2,"test SyncSlider, position is correct");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr4","test SyncSlider, correct chromosome");
        is($w->{count_covered},0,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{middle},2,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),0,"Correct number of data entries");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr4","test SyncSlider, correct chromosome");
        is($w->{count_covered},1,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test SyncSlider, correct number of snps in region");
        is($w->{start},1,"test SyncSlider, correct start position");
        is($w->{end},4,"test SyncSlider, correct end position");
        is($w->{middle},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),1,"Correct number of data entries");
        
        $w=$bpsl->nextWindow();
        not_exists($w,"correct end of file reached");
        $w=$bpsl->nextWindow();
        not_exists($w,"correct end of file reached");
        
        
        # weired characters
        $teststr=
        "chr1\t1\tN\t0:0:6:0:1:0\t0:0:0:7:3:0\t0:0:8:0:2:0\n".
        "chr1\t2\tN\t0:4:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:2\n".
        "chr1\t3\tN\t0:0:3:0:0:0\t0:0:0:8:0:0\t0:0:13:0:2:2\n".
        "chr1\t4\tN\t0:4:8:0:0:0\t-\t0:0:8:0:2:1\n";
        $bpsl=_getSyncSliderForString($teststr,3,1,2,4);
        #str, window, step, mincount, mincov, maxcov
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr1","test SyncSlider, correct chromosome");
        is($w->{count_covered},2,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is($w->{data}[0]{issnp},1,"test SyncSlider, snp correct");
        is($w->{data}[0]{ispuresnp},1,"test SyncSlider, snp correct");
        is($w->{data}[0]{iscov},1,"test SyncSlider, snp correct");
        is($w->{data}[1]{issnp},1,"test SyncSlider, snp correct");
        is($w->{data}[1]{ispuresnp},0,"test SyncSlider, snp correct");
        is($w->{data}[1]{iscov},1,"test SyncSlider, snp correct");
        is($w->{data}[2]{issnp},0,"test SyncSlider, snp correct");
        is($w->{data}[2]{ispuresnp},0,"test SyncSlider, snp correct");
        is($w->{data}[2]{iscov},0,"test SyncSlider, snp correct");
        
        is($w->{data}[0]{samples}[0]{eucov},6,"test SyncSlider, coverage is ok");
        is($w->{data}[0]{samples}[0]{totcov},7,"test SyncSlider, coverage is ok");
        is($w->{data}[0]{samples}[1]{eucov},7,"test SyncSlider, coverage is ok");
        is($w->{data}[0]{samples}[1]{totcov},10,"test SyncSlider, coverage is ok");
        
        is($w->{data}[1]{samples}[1]{eucov},8,"test SyncSlider, coverage is ok");
        is($w->{data}[1]{samples}[1]{totcov},8,"test SyncSlider, coverage is ok");
        is($w->{data}[1]{samples}[2]{eucov},8,"test SyncSlider, coverage is ok");
        is($w->{data}[1]{samples}[2]{totcov},10,"test SyncSlider, coverage is ok");
        is($w->{data}[1]{samples}[2]{del},2,"test SyncSlider, deletion count is ok");
        is($w->{data}[1]{samples}[1]{del},0,"test SyncSlider, deletion count is ok");
        is($w->{data}[1]{samples}[2]{A},0,"test SyncSlider, A count is ok");
        is($w->{data}[1]{samples}[2]{T},0,"test SyncSlider, T count is ok");
        is($w->{data}[1]{samples}[2]{C},8,"test SyncSlider, C count is ok");
        is($w->{data}[1]{samples}[2]{G},0,"test SyncSlider, G count is ok");
        is($w->{data}[1]{samples}[2]{N},0,"test SyncSlider, N count is ok");
        is($w->{data}[0]{samples}[2]{N},2,"test SyncSlider, N count is ok");
        
        
        $w=$bpsl->nextWindow();
        is($w->{data}[2]{pos},4,"test SyncSlider, position correct");
        is($w->{data}[2]{issnp},0,"test SyncSlider, snp correct");
        is($w->{data}[2]{ispuresnp},0,"test SyncSlider, snp correct");
        is($w->{data}[2]{iscov},0,"test SyncSlider, snp correct");
        is($w->{data}[2]{samples}[1]{del},0,"test SyncSlider, deletion count is ok");
        is($w->{data}[2]{samples}[1]{A},0,"test SyncSlider, A count is ok");
        is($w->{data}[2]{samples}[1]{T},0,"test SyncSlider, T count is ok");
        is($w->{data}[2]{samples}[1]{C},0,"test SyncSlider, C count is ok");
        is($w->{data}[2]{samples}[1]{G},0,"test SyncSlider, G count is ok");
        is($w->{data}[2]{samples}[1]{N},0,"test SyncSlider, N count is ok");
        is($w->{data}[2]{samples}[2]{del},1,"test SyncSlider, deletion count is ok");
        is($w->{data}[2]{samples}[2]{A},0,"test SyncSlider, A count is ok");
        is($w->{data}[2]{samples}[2]{T},0,"test SyncSlider, T count is ok");
        is($w->{data}[2]{samples}[2]{C},8,"test SyncSlider, C count is ok");
        is($w->{data}[2]{samples}[2]{G},0,"test SyncSlider, G count is ok");
        is($w->{data}[2]{samples}[2]{N},2,"test SyncSlider, N count is ok");
        
        
        # SNP calling as sum of its parts
        $teststr=
        "chr1\t1\tN\t3:1:0:0:0:0\t3:1:0:0:0:0\t4:0:0:0:0:0\n".
        "chr1\t2\tN\t4:0:0:0:0:0\t2:2:0:0:0:0\t4:0:0:0:0:0\n".
        "chr1\t3\tN\t3:1:0:0:0:0\t4:0:0:0:0:0\t4:0:0:0:0:0\n";
        $bpsl=_getSyncSliderForString($teststr,3,1,2,4);
        $w=$bpsl->nextWindow();
        
        is($w->{data}[0]{pos},1,"test SyncSlider, position correct");
        is($w->{data}[0]{issnp},1,"test SyncSlider, snp correct");
        is($w->{data}[0]{ispuresnp},1,"test SyncSlider, snp correct");
        is($w->{data}[0]{iscov},1,"test SyncSlider, coverage correct");
        
        is($w->{data}[1]{pos},2,"test SyncSlider, position correct");
        is($w->{data}[1]{issnp},1,"test SyncSlider, snp correct");
        is($w->{data}[1]{ispuresnp},1,"test SyncSlider, snp correct");
        is($w->{data}[1]{iscov},1,"test SyncSlider, coverage correct");
        
        is($w->{data}[2]{pos},3,"test SyncSlider, position correct");
        is($w->{data}[2]{issnp},0,"test SyncSlider, snp correct");
        is($w->{data}[2]{ispuresnp},0,"test SyncSlider, snp correct");
        is($w->{data}[2]{iscov},1,"test SyncSlider, coverage correct");
    }
    
}
1;