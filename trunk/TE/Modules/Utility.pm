{
    package Utility;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";

    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(load_te_annotation get_fasta_reader);
    
    sub load_te_annotation
    {
        my $file=shift;
        
        open my $ifh, "<", $file or die "could not open input file";
        
        my $h={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            
            my($id,$ann)=split /\t/,$l;
            
            $h->{$id}=$ann;
        }
        close $ifh;
        return $h;
    }
    
    sub get_fasta_reader
    {
        my $file=shift;
        open my $ifh,"<",$file or die "Could not open input file";
        my $lastheader=\"";
        
        return sub
        {
            my $header="";
            my $sequence="";
            while(my $l=<$ifh>)
            {
                chomp $l;
                if($l=~m/^>/)
                {
                    my($newheader)=$l=~m/^>(.+)/;
                    if($$lastheader)
                    {
                        $header=$$lastheader;
                        $lastheader=\$newheader;
                        return {
                            head=>$header,
                            seq=>$sequence
                        };
                    }
                    else
                    {
                        $lastheader=\$newheader;
                    }
                    
                }
                else
                {
                    $sequence.=$l;
                }
            }
            if($sequence)
            {
                return
                {
                  head=>$$lastheader,
                  seq=>$sequence
                };
            }
            else
            {
                return undef;
            }
            
        }
        
    }
    
}
1;