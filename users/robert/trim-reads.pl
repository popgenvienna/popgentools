{
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use IO::Compress::Gzip;
#use Pileup;
    my $minLength=0;

    my $input="";
    my $output="";

    
    GetOptions(
        "input=s"               =>\$input,
        "length=i"              =>\$minLength,
        "output=s"              =>\$output

    ) or pod2usage(-msg=>"Options not valid $!",-verbose=>1);
    
    # dive into the alternative branches (if requested)


    pod2usage(-msg=>"Length must be larger than 0",-verbose=>1) if $minLength<1 ;  # min length has to be 1 or larger
    pod2usage(-msg=>"At least one input file has to be provided", -verbose=>1) unless -e $input;
    pod2usage(-msg=>"An output file has to be provided", -verbose=>1) unless $output;
    


$output=$output . ".gz";
my  $ofh = new IO::Compress::Gzip $output or die "Could not open gzipped output file $output $!"; 

            
my $fastqr=FastqReader->new($input);
while(1)
        {
        my($firstheader,$nucleotide,$secondheader,$quality)=$fastqr->nextRead();
        last unless $firstheader;
        $nucleotide=substr($nucleotide,0,$minLength);
        $quality=substr($quality,0,$minLength);
        printFastq($ofh,$firstheader,$nucleotide,$secondheader,$quality);
        }
    

    sub printFastq
        {
            my $ofh=shift;
            my $header1=shift;
            my $nuc=shift;
            my $header2=shift;
            my $qual=shift;
            
            print $ofh $header1."\n";
            print $ofh $nuc."\n";
            print $ofh $header2."\n";
            print $ofh $qual."\n";
            
        }
        exit;
}


    
{
    package FastqReader;
    use IO::Uncompress::Gunzip;
    
    sub new
    {
        my $class=shift;
        my $file=shift;
        my $ofh = undef;
        if ($file=~/\.gz$/i) {
             $ofh = new IO::Uncompress::Gunzip $file or die "Could not open file gzipped file $file  $!";
	}
	else {
             open $ofh, "<", $file  or die "Could not open file handle, $!";
	}
    
        return bless {
            file=>$file,
            fh=>$ofh,
            buffer=>[]
        },__PACKAGE__;
    }

    
    sub nextRead
    {
        
        my $self=shift;

        my $firstheader=$self->nextLine;
        my $nucleotide=$self->nextLine;
        my $secondheader=$self->nextLine;
        my $quality=$self->nextLine;
        
        return undef unless $quality;        
        chomp $firstheader;
        chomp $nucleotide;
        chomp $secondheader;
        chomp $quality;
        
        die "first header must start with an @; line: $firstheader" unless $firstheader =~ m/^[@]/;
        die "second header must start with an +; line: $secondheader" unless $secondheader=~m/^[+]/;
        
        return ($firstheader,$nucleotide,$secondheader,$quality);
    }
    
    sub nextLine
    {
        my $self=shift;
        my $fh=$self->{fh};
        
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub bufferLine
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }


}




=head1 NAME

trim-fastq.pl - Perform quality filtering of fastq files

=head1 SYNOPSIS

 # Minimum argument call; single read mode
 trim-fastq.pl --input1 input.fastq --output output.fastq
 
# Minimum argument call; single read mode
 trim-fastq.pl --input1 input_1.fastq --input2 input_2.fastq --output output_prefix

=head1 OPTIONS

=over 4

=item B<--input1>

The input file, or the input file of the first read, in fastq format. Mandatory parameter

=item B<--output>

The output file of the first read. Will be in fastq. Mandatory parameter

=item B<--length>

The minimum length of the read after trimming; default=40

=item B<--quit>

suppress output to stdout (console)

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS

The script removes 'N' - characters at the beginning and the end of the provided reads. If any remaining 'N' characters are found the read is discarded.
Quality removal is done using a modified Mott-algorithm;
For each base a score is calculated: score_base = quality_base - threshold

While scanning along the read a running sum of this score is calculatet; If the score drops below zero the score is set to zero;
The highest scoring region of the read is finally reported;

=head2 INPUT

Input must be a fastq file. No line breaks in the nucleotide or quality sequence are allowed. No empty lines in the file are allowed.
The fastq file produced by the Illumina pipeline can be directyl use. Example input:

    @read1
    NNNNNAAAAAAAAAATTTTTTTTTTAAAAAAAAAANNNNNNNNNN
    +read1
    UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
    @read2
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT
    +read2
    UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUSSSSSSSSSS



=head2 OUTPUT

the output file will be in fastq. The quality sequence will be provided in the same format as in the input file (apart from trimming);
    
    @read1
    AAAAAAAAAATTTTTTTTTTAAAAAAAAAA
    +read1
    UUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
    @read2
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    +read2
    UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU

=cut