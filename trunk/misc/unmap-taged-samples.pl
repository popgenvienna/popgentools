use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";


#my @bc=qw/GACTGGACT TCAGCAGT CTGATCT AGTCAT/;

my $fastqfile1;
my $fastqfile2;
my $rawtagstring;
my $bothreads_taged=0;
my $maxmismatches=0;
my $output="";
my $help=0;

GetOptions(
    "read1=s"               =>\$fastqfile1,
    "read2=s"               =>\$fastqfile2,
    "tags=s"                =>\$rawtagstring,
    "output=s"              =>\$output,
    "max-mismatches=i"       =>\$maxmismatches,
    "tag-in-both-reads"     =>\$bothreads_taged,   
    "help"                  =>\$help

) or die "Wrong parameters";
pod2usage(-verbose=>2) if $help;
pod2usage(-verbose=>1,-message=>"You have to specify the first fastq file") unless -e $fastqfile1;
pod2usage(-verbose=>1,-message=>"You have to specify the second fastq file") unless -e $fastqfile2;
pod2usage(-verbose=>1,-message=>"You have to provide a output file") unless $output;
pod2usage(-verbose=>1,-message=>"You have to provide some tags") unless $rawtagstring;


my $tags=Utility::load_tags($rawtagstring,$output);
Utility::check_tags($tags,$maxmismatches);

my $fr1 = Utility::get_fastq_reader($fastqfile1);
my $fr2 = Utility::get_fastq_reader($fastqfile2);

my $istagmatcher=undef;
if($maxmismatches)
{
    $istagmatcher=\&Utility::istagmatch_mismatches;
}
else
{
    $istagmatcher=\&Utility::istagmatch_perfect;
}

while(1)
{
    my $s1=$fr1->();
    my $s2=$fr2->();
    last unless $s1;
    
    # identifiy tags
    # my($nr1,$nr2)=("","");
    my $tr1=undef;
    my $tr2=undef;

    foreach my $tag (@$tags)
    {
        $tr1 = $tag if $istagmatcher->($s1->{seq},$tag,$maxmismatches);
        last if $tr1;
    }
    
    foreach my $tag(@$tags)
    {
        $tr1=$tag if $istagmatcher->($s2->{seq},$tag,$maxmismatches);
        last if $tr2;
    }
    

    # store according to tag info
    if($tr1 and $tr2)
    {
        # are both tags identified?
        if($tr1->{tag} eq $tr2->{tag})
        {
            # are they identical?

            Utility::write_sequences($s1,$s2,$tr1);
        }
    }
    elsif($tr1 or $tr2)
    {
        unless($bothreads_taged)
        {
            # is at least one tag identified: if so used it as well
            my $tr =$tr1?$tr1:$tr2;
            Utility::write_sequences($s1,$s2,$tr);
        }
    }
}



exit;

{
    package Utility;
    use strict;
    use warnings;
    
    sub istagmatch_perfect
    {
        my $sequence=shift;
        my $tag=shift;
        my $mismatches=shift;
        die "Only zero mismatches allowed for perfect matching" if $mismatches;
        my $tagseq=$tag->{tag};
        return 1 if $sequence=~m/^$tagseq/;
        return 0;
    }
    
    sub istagmatch_mismatches
    {
        my $sequence=shift;
        my $tag=shift;
        my $mismatches=shift;
        my @seq=split //,$sequence;
        my @tagar=@{$tag->{tagar}};
        
        my $dist=0;
        for(my $i=0; $i<@tagar; $i++)
        {
            $dist++ if $seq[$i] ne $tagar[$i];
        }
        
        return 1 if $dist< $mismatches;
        return 0;
    }
    
    sub check_tags
    {
        my $tags=shift;
        my $maxmismatches=shift;
        my $tagcount=@$tags;        
        die "At least two tags need to be provided" unless $tagcount>1;


        my $mineditdist=undef;
        foreach(my $i=0; $i<$tagcount; $i++)
        {
            foreach(my $k=$i+1; $k<$tagcount; $k++)
            {
                my $tag1=$tags->[$i];   #tag, leng, f1, f2, fh1, fh2, tagar
                my $tag2=$tags->[$k];
                
                my @tar1=@{$tag1->{tagar}};
                my @tar2=@{$tag2->{tagar}};
                
                my $editdist=0;            
                my $minleng= scalar(@tar1) <= scalar(@tar2) ? scalar(@tar1) : scalar(@tar2);
                for(my $j=0; $j<$minleng; $j++)
                {
                    if($tar1[$j] ne $tar2[$j])
                    {
                        $editdist++
                    }
                }
                $mineditdist=$editdist if(not defined($mineditdist));
                $mineditdist=$editdist if($editdist<$mineditdist);
            }
        }
        my $tocomp=$mineditdist/2;
        if($maxmismatches<$tocomp)
        {
            print "OK: number of allomed mismatches ($maxmismatches) is smaller than half the minimum edit distance between any pair of sequences ($mineditdist/2)\n"
        }
        else
        {
            die "Can not resolve tags with the given number of allowed mismatches\nHalf of the minimum edit distance between the tags ($mineditdist/2) needs to be larger than the number of allowed mismatches ($maxmismatches)";
        }

    }
    
    sub load_tags
    {
        my $rawtag=shift;
        my $output=shift;
        my @tags=split /\s+/,$rawtag;
        push @tags, $rawtag unless @tags;
        
        
        my $tc=@tags;
        for my $i(0..($tc-1))
        {
            for my $k(($i+1)..($tc-1))
            {
                my $t1=$tags[$i];
                my $t2=$tags[$k];
                die"tag $t1 and $t2 have the same starting sequence, invalid!" if($t1=~/^$t2/ or $t2=~/^$t1/)
            }
        }
        
        my $tagar=[];
        foreach my $t(@tags)
        {
            my $op1=$output."_".$t."_1";
            my $op2=$output."_".$t."_2";
            open my $fh1,">", $op1 or die "Could not open file handle";
            open my $fh2,">", $op2 or die "Could not open file handle";
            push @$tagar,{
                tag=>$t,
                leng=>length($t),
                tagar=>[split //,$t],
                f1=>$op1,
                f2=>$op2,
                fh1=>$fh1,
                fh2=>$fh2
            };
        }
        return $tagar; #tag, leng, f1, f2, fh1, fh2, tagar
    }
    
    sub write_sequences
    {
        my $s1=shift;
        my $s2=shift;
        my $tag=shift;
        
        my $leng=$tag->{leng};
        $s1->{seq}=~s/^.{$leng}//;
        $s1->{qual}=~s/^.{$leng}//;
        $s2->{seq}=~s/^.{$leng}//;
        $s2->{qual}=~s/^.{$leng}//;
        
        my $fh1=$tag->{fh1};
        my $fh2=$tag->{fh2};
        
        print $fh1 $s1->{h1}."\n";
        print $fh1 $s1->{seq}."\n";
        print $fh1 $s1->{h2}."\n";
        print $fh1 $s1->{qual}."\n";

        print $fh2 $s2->{h1}."\n";
        print $fh2 $s2->{seq}."\n";
        print $fh2 $s2->{h2}."\n";
        print $fh2 $s2->{qual}."\n";        
        
    }
    
    sub get_fastq_reader
    {
        my $file=shift;
        open my $ifh, "<",$file or die "Could not open input file";
        
        return sub
        {
            my $h1=<$ifh>;
            my $seq=<$ifh>;
            my $h2=<$ifh>;
            my $qual=<$ifh>;
            return undef unless $h1;
            chomp $seq; chomp $h1; chomp $qual; chomp $h2;
            return
            {
                h1=>$h1,
                seq=>$seq,
                h2=>$h2,
                qual=>$qual
            };
        }
    }
}

    #"read1=s"               =>\$fastqfile1,
    #"read2=s"               =>\$fastqfile2,
    #"tags=s"                =>\$rawtagstring,
    #"output=s"              =>\$output,
    #"help"                  =>\$help

=head1 NAME

unmap-taged-samples.pl 

=head1 SYNOPSIS

 unmap-taged-samples --read1 sequence_1.txt --read_2 sequence_2.txt --tags "ATCG GATCA TCAGTA GGTGATA" --output unmaped --max-mismatches 2

=cut