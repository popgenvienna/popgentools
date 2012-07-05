{
package TESiteUtility;
use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use TESite;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(load_te_sites write_te_sites);

sub load_te_sites
{
    my $te_insert_file=shift;
    open my $ifh, "<", $te_insert_file or die "Could not open input file";
    my $te_inserts=[];
    while(my $l=<$ifh>)
    {
        chomp $l;
        next if $l=~/^\s*$/;
        my $tei=TESite->new(split /\t/,$l);
        # TESite -> # chr, insdir, teid, support, start, end, comment
        die "No TE site" unless $tei; 
        push @$te_inserts,$tei;
    }
    return $te_inserts;
}


sub write_te_sites
{
    my $file=shift;
    my $tesites=shift;
    open my $ofh ,">",$file or die "Could not open output file";
    foreach my $ts (@$tesites)
    {
        my $form=$ts->tostring();
        print $ofh, $form."\n";
    }
}


}
1;
