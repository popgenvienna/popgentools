package Module;
use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @Export=qw();
our @EXPORT_OK = qw();

sub load_divergence
{
    my $file=shift;
    my $annotation=shift;
        #YHet    291809  T       A       T
        #YHet    291810  T       T       T
        #YHet    291811  T       A       A
        #YHet    291812  T       C       T
    
    open my $ifh, "<", $annotation or die "Could not open input file";
    my $l=<$ifh>;
    close $ifh;
    chomp $l;
    my @ar=split /\t/,$l;
    my $columncount=@ar;
    if($columncount==5)
    {
        return _load_divergence_withoutgroup($file,$annotation);
    }
    elsif($columncount==4)
    {
        return _load_pairwise_divergence($file,$annotation);
    }
    else
    {
        die "Divergence file does not meet the requirements: $l\n";
    }
    

    
}

sub _load_divergence_withoutgroup
{
    my $file=shift;
    my $annotation=shift;
        #YHet    291809  T       A       T
        #YHet    291810  T       T       T
        #YHet    291811  T       A       A
        #YHet    291812  T       C       T
    
    my $toret={};
    open my $ifh, "<", $file or die "Could not open input file";
    while(my $l=<$ifh>)
    {
        chomp $l;
        my($chr,$pos,$ref,$comp,$outgroup)=split /\t/,$l;
        # kick out everything which is not annotated
        next unless(exists($annotation->{$chr}{$pos}));
        next if($annotation->{$chr}{$pos}==7);
        my $diverged=0;
        if($ref ne $comp and $ref ne $outgroup and $comp ne $outgroup)
        {
            # A T C
            $diverged=0;
        }
        else
        {
            if($ref eq $comp or $ref eq $outgroup)
            {
                # A A T
                # A A A
                # A T A
                $diverged=$ref.$ref;
            }
            elsif($ref ne $comp and $comp eq $outgroup)
            {
                # T A A
                $diverged=$ref.$comp
            }
            else
            {
                die "unconsidered case of divergence $l";
            }
        }
        $toret->{$chr}{$pos}=$diverged;
    }
    return $toret;
    
}

sub _load_pairwise_divergence
{
    my $file=shift;
    my $annotation=shift;
        #YHet    291809  T       A
        #YHet    291810  T       T
        #YHet    291811  T       A
        #YHet    291812  T       C 
    
    my $toret={};
    open my $ifh, "<", $file or die "Could not open input file";
    while(my $l=<$ifh>)
    {
        chomp $l;
        my($chr,$pos,$ref,$comp)=split /\t/,$l;
        # kick out everything which is not annotated
        next unless(exists($annotation->{$chr}{$pos}));
        next if($annotation->{$chr}{$pos}==7);
        my $diverged;
        if($ref ne $comp)
        {
            # A T
            $diverged=$ref.$comp;
        }
        elsif($ref eq $comp)
        {
            $diverged=$ref.$ref;
        }
        $toret->{$chr}{$pos}=$diverged
    }
    return $toret;
}





1;
