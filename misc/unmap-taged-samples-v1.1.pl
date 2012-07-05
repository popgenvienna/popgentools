use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw/$RealBin/;
#use lib "$RealBin/../Modules";

# Note: This script search for first 5 bp of tags but will trim the whole barcode length from 5' end.
# Modified date: 9th Jan, 2012

#my @bc=qw/GACTGGACT TCAGCAGT CTGATCT AGTCAT/;

my $fastqfile1;
my $fastqfile2;
my $rawtagstring;
my $bothreads_taged=0;
my $output="";
my $help=0;

GetOptions(
    "read1=s"               =>\$fastqfile1,
    "read2=s"               =>\$fastqfile2,
    "tag-file=s"                =>\$rawtagstring, # tag tab delimeted file
    "output=s"              =>\$output,
    "tag-in-both-reads"     =>\$bothreads_taged,   
    "help"                  =>\$help

) or die "Wrong parameters";
pod2usage(-verbose=>2) if $help;
pod2usage(-verbose=>1,-message=>"You have to specify the first fastq file") unless -e $fastqfile1;
pod2usage(-verbose=>1,-message=>"You have to specify the second fastq file") unless -e $fastqfile2;
pod2usage(-verbose=>1,-message=>"You have to provide a output file") unless $output;
pod2usage(-verbose=>1,-message=>"You have to provide some tags") unless $rawtagstring;

my $tags = {};
my $partial_tag_hash = {};
($tags,$partial_tag_hash)=Utility::load_tags($rawtagstring,$output);
my $fr1 = Utility::get_fastq_reader($fastqfile1);
my $fr2 = Utility::get_fastq_reader($fastqfile2);

foreach my $key (keys %$tags) {
    #print "$key=>$tags->{$key}->{leng}\n";
}
my @bc=keys(%$tags);

die "The third column barcodes are not unique!!\n" unless(scalar(keys %$partial_tag_hash) == scalar(keys %$tags));

#exit();
while(1)
{
    my $s1=$fr1->();
    my $s2=$fr2->();
    last unless $s1;
    
    # identifiy tags
    my($nr1,$nr2)=("","");
    foreach my $bc (keys %$partial_tag_hash)
    {
        $nr1=$bc if $s1->{seq}=~/^$bc/;
        $nr2=$bc if $s2->{seq}=~/^$bc/;
    }
    

    # store according to tag info
    if($nr1 and $nr2)
    {
        # are both tags identified?
        if($nr1 eq $nr2)
        {
            # are they identical?
            #my $tag=$nr1;
            my $tag=$partial_tag_hash->{$nr1};
            #print "hi1,,$tag,,$tags->{$tag}->{leng}\n";
            Utility::write_sequences($s1,$s2,$tags->{$tag});
        }
    }
    elsif($nr1 or $nr2)
    {
        unless($bothreads_taged)
        {
            
            # is at least one tag identified: if so used it as well
            my $tag =$nr1?$nr1:$nr2;
            $tag=$partial_tag_hash->{$tag};
            #print "hi2,,$tag,,$tags->{$tag}->{leng}\n";
            Utility::write_sequences($s1,$s2,$tags->{$tag});
        }
    }
}



exit;

{
    package Utility;
    use strict;
    use warnings;
    
    sub load_tags
    {
        my $rawtag=shift;
        my $output=shift;
        

        my @tags = ();
        my @tags1 = ();

        open my $tagfh,"<", $rawtag or die "Could not open file handle";
        while(<$tagfh>) {
            chomp;
	    my $line = $_;

	    # discard blank line
	    if ($line =~ m/^\s*$/g) {
		next;
	    }
	    # discard comment line
	    elsif($line =~ m/^\s*#/g) {
		next;
	    }
	    else {
                my ($id,$barcode_actual,$barcode_initial) = ("","","");
                my @line = split("\t",$line);
                
                $id = shift(@line);
                $barcode_actual = shift(@line);
                $barcode_initial = shift(@line);
                push(@tags,$barcode_actual);
                push(@tags1,"$barcode_actual#$barcode_initial");
            }
        }
        
        close $tagfh;
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
        
        my $taghash={};
        my $partial_tag_hash={};
        foreach my $t(@tags1)
        {
            my ($barcode_actual,$barcode_initial) = ("","");
            ($barcode_actual,$barcode_initial) = split("#",$t);
            
            $partial_tag_hash->{$barcode_initial} = $barcode_actual;
            
            $t = $barcode_actual;
            my $op1=$output."_".$t."_1";
            my $op2=$output."_".$t."_2";
            open my $fh1,">", $op1 or die "Could not open file handle";
            open my $fh2,">", $op2 or die "Could not open file handle";
            $taghash->{$t}={
                tag=>$t,
                leng=>length($t),
                f1=>$op1,
                f2=>$op2,
                fh1=>$fh1,
                fh2=>$fh2
            };
        }
        return ($taghash,$partial_tag_hash); #tag, leng, f1, f2, fh1, fh2
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


=head1 NAME

unmap-taged-samples-v1.1.pl

=head1 SYNOPSIS

 perl unmap-taged-samples-v1.1.pl -read1 sequence_1.txt --read_2 sequence_2.txt --tag-file tag.txt --output unmaped

=head1 OPTIONS

=over 4

=item B<--read1>

The read1 fastq file. Mandatory parameter

=item B<--read2>

The read2 fastq file. Mandatory parameter

=item B<--tag-file>

The barcode tab delemeited file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

A tab delemited tag file; example:

 PEMx5	    GACTGGACT	GACTG
 PEMx6	    TCAGCAGT	TCAGC
 PEMx7	    CTGATCT	CTGAT
 PEMx8	    AGTCAT	AGTCA
 PEMx9	    CGTACTGAT	CGTAC
 PEMx10	    ATGCTACT	ATGCT
 PEMx11	    GCATGCT	GCATG
 PEMx12	    TACGAT	TACGA
 PEMx13	    GAGACGTCT	GAGAC
 PEMx14	    AGAACGAT	AGAAC
 PEMx15	    CAGTAAT	CAGTA
 PEMx16	    ATACGT	ATACG

 col 1: sample id
 col 2: full barcode to be clipped
 col 3: part of full barcode to be searched for clipping the full barcode

=head1 AUTHORS

Ram vinay pandey
Robert Kofler


=cut

