use strict;
use warnings;

my $filename=shift;
my $outputfile=shift;
open my $ifh, "<", $filename or die "could not open input file\n";

die "no output file" unless $outputfile;

# __fig_outline_method__
# __sup_tefreqspec__
# __tab_posselcandidates__

# file slurping
undef $/;
my $file =<$ifh>;


my $resolvh=get_substitute_hash($file);

while($file=~m/(__.*?__)/)
{
    my $m=$1;
    my $key=$m;
    $key=~s/\s//;
    
    my $value=$resolvh->{$key};
    die "no value for $key, $m" unless $value;
    $file=~s/$m/$value/;
}


open my $ofh,">",$outputfile or die "could not open output file";
print $ofh $file;
close $ofh;

while(my($key,$value)=each(%$resolvh))
{
    print "$key\t$value\n";
}

exit;

sub get_substitute_hash
{
    my $file=shift;
    my $tabl=[];
    my $figl=[];
    my $supl=[];
    
    while($file=~m/(__.*?__)/g)
    {
        my $m=$1;
        $m=~s/\s+//g;
        if($m=~/^__fig/)
        {
            push @$figl,$m;
        }
        elsif($m=~/^__tab/)
        {
            push @$tabl,$m;
        }
        elsif($m=~/^__sup/)
        {
            push @$supl,$m;
        }
        else
        {
            my $before =$`;
            my $after=$';
            $before=substr($before,length($before)-20,20);
            $after=substr($after,0,20);
            warn "Could not resolve $m;\n $before:$m:$after\n";
        }
    }
    
    my $resolvh={};
    my $counter=1;
    foreach my $t (@$tabl)
    {
        next if(exists($resolvh->{$t}));
        $resolvh->{$t}=$counter;
        $counter++;
    }
    
    $counter=1;
    foreach my $f (@$figl)
    {
        next if(exists($resolvh->{$f}));
        $resolvh->{$f}=$counter;
        $counter++;
    }
    
    $counter=1;
    foreach my $s (@$supl)
    {
        next if(exists($resolvh->{$s}));
        $resolvh->{$s}=$counter;
        $counter++;
    }
    return $resolvh;
    
}