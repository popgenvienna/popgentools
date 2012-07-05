use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use Pileup;

my $encoder=get_quality_encoder("illumina");

my $file=shift;
my $fastqr=FastqReader->new($file);


my($lengsum,$qualsum)=(0,0);
my($readcount)=(0,0);

    while(1)
    {
        my($firstheader,$nucleotide,$secondheader,$quality)=$fastqr->nextRead();
        last unless $firstheader;
        
        $readcount++;
        $lengsum+=length($nucleotide);
        
        my @quals=split//,$quality;
        foreach my $q(@quals)
        {
            my $q=$encoder->($q);
            $qualsum+=$q;
        }
     
    }
    
    my $avleng=$lengsum/$readcount;
    my $avqual=$qualsum/$lengsum;
    print "Reads:\t$readcount\n";
    print "Bases:\t$lengsum\n";
    print "Average length\t$avleng\n";
    print "Average quality\t$avqual\n";
    exit;
    
    
{
    package FastqReader;
    
    sub new
    {
        my $class=shift;
        my $file=shift;
        open my $ofh,"<$file" or die "Could not open file handle";
        
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