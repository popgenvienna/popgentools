{
package TESite;
use strict;
use warnings;

# chr, insdir, teid, support, start, end, comment

sub new {
    my $class = shift;
    my $chr=shift;
    my $insdir=shift;
    my $teid=shift;
    my $support=shift;
    my $start=shift;
    my $end=shift;
    my $comment=shift;
    my $self = bless {
                      chr       =>$chr,
                      insdir    =>$insdir,
                      teid      =>$teid,
                      support   =>$support,
                      start     =>$start,
                      end       =>$end,
                      comment   =>$comment
                      }, __PACKAGE__;
    return $self;
}


sub key
{
    my $self=shift;
    return "$self->{chr}:$self->{insdir}:$self->{teid}:$self->{start}:$self->{end}";
}

sub tostring
{
    my $self=shift;
    my @ar=($self->{chr}, $self->{insdir}, $self->{teid}, $self->{support}, $self->{start}, $self->{end}, $self->{comment});
    return join("\t",@ar);
    
    
}

}


1;





