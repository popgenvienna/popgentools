{
    package SynchronizeUtility;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT =qw(format_parsed_pileup);
    
    
    sub format_parsed_pileup
    {
        my $pp=shift;
        return "$pp->{A}:$pp->{T}:$pp->{C}:$pp->{G}:$pp->{N}:$pp->{del}" if $pp;
        return "-";
    }
    


}

1;