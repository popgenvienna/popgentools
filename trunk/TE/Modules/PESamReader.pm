{
    package PESamReader;
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";

    
    
    sub new
    {
        my $class=shift;
        my $file=shift;

        
        open my $fh,"<",$file or die "Could not open file handle";
        
        
        my $tobless={
            file    =>$file,
            pair    =>{},
            fh      =>$fh,
            buffer  =>[]
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
    

    
    
    sub next
    {
        my $self=shift;

        my $ph=$self->{pair};
        
        while(1)
        {
            my $line=$self->_nextline();
            last unless $line;
            chomp $line;
            my $id=$self->_getID($line);
            
            if(exists($ph->{$id}))
                {
                    my $lastline=$ph->{$id};
                    my $novelline=$line;
                    delete($ph->{$id}); # delete the hash entry -> avoid memory overkill
                    return [$line,$lastline];
                }
                else
                {
                    $ph->{$id}=$line;
                }
            
        }
        return undef;
    }
    
    sub _getID{
        my $self=shift;
        my $line=shift;
        my @a= split /\t/, $line;
        my $id=$a[0];
        return $id;
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
    
    sub close
    {
        my $self=shift;
        my $fh=$self->{fh};
        close $fh;
    }
}

1;