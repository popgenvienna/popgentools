#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use List::Util qw[min max];
use TEHierarchy;
use ParseSam;

my $input;
my $tehierfile;
my $tehiertargetlevel;
my $mincount=2;
my $insertdistance=80;
my $insertdistStddev=35;
my $readlength=75;
my $narrowrange=75;
my $output="";
my $minMapquality=20;
my $help=0;
my $test=0;

#--input Volumes/Volume_3/analysis/te/ky_es_po/map/Es_Ky_Po.sam --te-hierarchy-file /Volumes/Volume_3/analysis/te/ky_es_po/ref/fb-hierarchy --te-hierarchy-level familiy --narrow-range 100 --min-count 3 --min-map-qual 15 --output /Volumes/Volume_3/analysis/te/ky_es_po/te-inserts/raw/EsKyPo.sites
# --input /Volumes/Volume_3/analysis/te/ky_es_po/map/EsKyPo_euchrom.sam --te-hierarchy-file /Volumes/Volume_3/analysis/te/ky_es_po/ref/fb-hierarchy --te-hierarchy-level familiy --narrow-range 100 --min-count 3 --min-map-qual 15 --output /Volumes/Volume_3/analysis/te/ky_es_po/te-inserts/raw/EsKyPo.sites

GetOptions(
    "input=s"	            =>\$input,
    "output=s"              =>\$output,
    "min-count=i"           =>\$mincount,
    "narrow-range=i"        =>\$narrowrange,
    "min-map-qual=i"        =>\$minMapquality,
    "te-hierarchy-file=s"   =>\$tehierfile,
    "te-hierarchy-level=s"  =>\$tehiertargetlevel,
    "test"                  =>\$test,
    "help"	            =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help; 
pod2usage(-verbose=>1,-msg=>"Input file does not exist") unless -e $input;
pod2usage(-verbose=>1,-msg=>"No output file provided") unless $output;
pod2usage(-verbose=>1,-msg=>"TE hierarchy file does not exist") unless -e $tehierfile;
pod2usage(-verbose=>1,-msg=>"No target level for TE hierarchy provided") unless $tehiertargetlevel;
 
my $teh_resolver=get_te_hierarchy_resolver($tehierfile,$tehiertargetlevel);

my $cp=Utility::get_chunk_parser(($readlength+$insertdistance+2*$insertdistStddev),$narrowrange,$mincount, $teh_resolver);
# 75 + 80 + 70 = 225

open my $ofh, ">", $output or die "Could not open output file";

my $samParser=get_te_samparser($teh_resolver);
my $scr=SamChunkReader->new($input,($readlength+$insertdistance+5*$insertdistStddev),$samParser,$minMapquality);

while(my $sc=$scr->nextChunk())
{
    next unless @$sc;
    my ($chr,$start,$end,$count)=($sc->[0]{chr},$sc->[0]{start},$sc->[-1]{end},scalar(@$sc));
    print "Processing Chunk of $count TE-insertion reads; $chr:$start-$end; ";
    
    # chr,strand, teid, count, nstart, nend, wstart, wend, posstring, sam
    my $annotes=$cp->($sc);
    foreach my $an (@$annotes)
    {
        print $ofh "$an->{chr}\t$an->{strand}\t$an->{teid}\t$an->{count}\t$an->{nstart}\t$an->{nend}\t$an->{wstart}:$an->{wend}=>$an->{posstring}\n"
    }
    my $counti=@$annotes;
    print "Identified $counti TE-insertion sites\n";


}
exit;

{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    sub get_chunk_parser
    {
        # get a chunk from the chromosome which contains varien TE insertions and resolve and annotate them
        my $wide_range=shift;  # all TEs of one kind which belong together (to the same insertion)
        my $narrow_range=shift; # the range which will be considered for polymorphism detection. Will be more narrow than the wide range
        my $mincount=shift;
        my $te_translator=shift;
        die "wide range must be larger than the narrow range" if $wide_range < $narrow_range;
        
        return sub
        {
            my $chunk=shift;
            my $teh={};
            foreach my $te (@$chunk)
            {
                # readid, flag, chr, start, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
                # end, end_s, start_s, te_ins_read, te_ins_direction
                my $str=$te->{te_ins_direction};
                my $chrmate=$te->{chrmate};
                my $teid=$te_translator->($chrmate);
                $teh->{$str}{$teid}=[] unless exists($teh->{$str}{$teid});
                push @{$teh->{$str}{$teid}},$te;
            }
            
            my $ranges=[];
            foreach my $strand (("F","R"))
            {
                while(my($te,$insertions)=each(%{$teh->{$strand}}))
                {
                        my $tes=_parseRange($insertions,$strand,$wide_range, $narrow_range, $mincount,$te);
                        push @$ranges,@$tes;
                }
            }
            return $ranges;    
        }
    }
    
    sub _parseRange
    {
        my $insertions=shift;
        my $strand=shift;
        my $wide_range=shift;
        my $narrow_range=shift;
        my $mincount=shift;
        my $teid=shift;
        my $feature=$strand eq "F"?"start_s": "end_s";
        
        @$insertions =sort {$a->{$feature}<=>$b->{$feature}} @$insertions;
    
        my $temp=[];
        my $lastpos=0;
        my $active=[];
        
        foreach my $te (@$insertions)
        {
            # the start position is crucial for the fwd insertions
            my $pos=$te->{$feature};
            $lastpos=$pos unless $lastpos;
            
            if($pos <$lastpos+$wide_range)
            {
                push @$active,$te;
            }
            else
            {
                push @$temp,$active;
                $active=[];
                push @$active,$te;
            }
            $lastpos=$pos;
        }
        push @$temp,$active if @$active;
        
        my $toret=[];
        foreach my $t(@$temp)
        {
            next unless @$t >=$mincount;
            $t = $strand eq "F" ? _annotateFwd($t,$narrow_range,$teid) : _annotateRev($t,$narrow_range,$teid);
            push @$toret,$t;
        }
        return $toret;
    }
    
    sub _annotateFwd
    {
        my $tes=shift;
        die "List of transposable elements must not be empty" unless @$tes;
        my $narrow_range=shift;
        my $teid=shift;
        
        my $chr=$tes->[0]{chr};
        my $minstart=min( map {$_->{start_s}} @$tes);
        my $maxstart=max(map{$_->{start_s}}@$tes);
        my $range=$maxstart-$minstart+1; 
        
        my $narrowstart=$maxstart-$narrow_range+1;
        $narrowstart=$minstart if $narrowstart<$minstart;

        # initialize the array,
        my $startposcount=[split //,"0" x $range];
        foreach my $te (@$tes)
        {
            my $relstart=$te->{start_s}-$minstart;
            $startposcount->[$relstart]++;
            
        }
        my $posstring=join(",",@$startposcount);
        
        # chr,strand, teid, count, nstart, nend, wstart, wend, posstring, sam
        my $entry=
        {
          chr=>$chr,
          strand=>"F",
          teid=>$teid,
          count=>scalar(@$tes),
          nstart=>$narrowstart,
          nend=>$maxstart,
          wstart=>$minstart,
          wend=>$maxstart,
          posstring=>$posstring,
          sam=>$tes
        };
        return $entry;
    }
    
    sub _annotateRev
    {
        my $tes=shift;
        die "List of transposable elements must not be empty" unless @$tes;
        my $narrow_range=shift;
        my $teid=shift;
        
        my $chr=$tes->[0]{chr};
        my $minend=min( map {$_->{end_s}} @$tes);
        my $maxend=max(map{$_->{end_s}}@$tes);
        my $range=$maxend-$minend+1; 
        
        my $narrowend=$minend+$narrow_range-1;
        $narrowend=$maxend if $narrowend>$maxend;

        # initialize the array,
        my $endposcount=[split //,"0" x $range];
        foreach my $te (@$tes)
        {
            my $relend=$te->{end_s}-$minend;
            $endposcount->[$relend]++;
        }
        my $posstring=join(",",@$endposcount);
        
        # chr, teid, count, nstart, nend, wstart, wend, posstring, sam
        my $entry=
        {
          chr=>$chr,
          strand=>"R",
          teid=>$teid,
          count=>scalar(@$tes),
          nstart=>$minend,
          nend=>$narrowend,
          wstart=>$minend,
          wend=>$maxend,
          posstring=>$posstring,
          sam=>$tes
        };
        return $entry;
    }  
    
}



{
    package SamChunkReader;
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use ParseSam;
    
    
    sub new
    {
        my $class=shift;
        my $file=shift;
        my $te_range=shift;  #suggest: readlength+ inner_dist + 3*stdev; max distance between start and next start
        my $sp=shift;
        my $minMapqual=shift;
        
        
        open my $fh,"<",$file or die "Could not open file handle";
        
        
        my $tobless={
            file    =>$file,
            fh      =>$fh,
            samparser=>$sp,
            buffer  =>[],
            minmapqual=>$minMapqual,
            range   =>$te_range     
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
    
    
    sub nextChunk
    {
        my $self=shift;
        # readid, flag, chr, start, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
        # end, end_s, te_ins_read, te_ins_direction
        my $sp = $self->{samparser};
        my $range = $self->{range};
        my $minmapqual=$self->{minmapqual};
        
        my $tereads=[];
        
        my $laststart=0;
        my $lastchr="";
        while(1)
        {
            my $line=$self->_nextline();
            last unless $line;
            chomp $line;
            my $s=$sp->($line);
            next if $s->{flag} &0x004;
            next if $s->{flag} &0x008;
            next unless ($s->{te_ins_read});
            next if ($s->{mq}<$minmapqual);
            $laststart  =$s->{start} unless $laststart;
            $laststart  =$s->{start} if($s->{start} <= $laststart+$range);
            $lastchr    =$s->{chr} unless $lastchr;
            
            if($s->{chr} ne $lastchr or $s->{start}> $laststart+$range)
            {
              $self->_bufferline($line);
              return $tereads;
            }
            else
            {
                push @$tereads, $s;   
            }
        }
        
        return $tereads if @$tereads;
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
}

#GetOptions(
#    "input=s"	            =>\$input,
#    "output=s"              =>\$output,
#    "min-count=i"           =>\$mincount,
#    "narrow-range=i"        =>\$narrowrange,
#    "min-map-qual=i"        =>\$minMapquality,
#    "te-hierarchy-file=s"   =>\$tehierfile,
#    "te-hierarchy-level=s"  =>\$tehiertargetlevel,
#    "test"                  =>\$test,
#    "help"	            =>\$help

=head1 NAME

perl identify-te-insertsites.pl - Identifies TE insertion sites (forward or reverse insertion) from a sam file

=head1 SYNOPSIS

 perl identify-te-insertsites.pl --input pe_maped_pool.sam --output te_insertionsites.txt --min-count 3 --narrow-range 100 --min-map-qual 15 --te-hierarchy-file fb-hierarchy.txt --te-hierarchy-level family

=head1 OPTIONS

=over 4

=item B<--input>

A sam file; Whole genome paired-end reads of a pooled population have to be mapped to a modified reference sequence consisting of repeat-masked reference chromosomes and the sequences of the transposons; Mandatory

=item B<--output>

The output file. Mandatory

=item B<--min-count>

the minimum number of PE-fragments that confirm the insertion of a TE of a certain family

=item B<--narrow-range>

the maximum length of the range which will be used for tallying the presence and absence fragments

=item B<--min-map-qual>

the minimum mapping quality; this will only apply to reads mapping to a reference contig.

=item B<--te-hierarchy-file>

a file containing a te hierarchy

=item B<--te-hierarchy-level>

the TE hierarchy level at which should be operated (eg.: family, order, class); depends on the categories present in the te-hierarchy-level

=item B<--help>

Display the help

=back

=head1 DETAILS

=cut