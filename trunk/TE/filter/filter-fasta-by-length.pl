#!/usr/bin/perl -w
use strict;
use warnings;

use lib '/Users/robertkofler/dev/PopGenTools/Modules';
use BasicUtility;

# USAGE: filter-fasta-by-length input-file.fasta minleng output-file 
my $input=shift;
my $output=shift;
my $minleng=shift;

my $fr=get_fasta_reader($input);
my $fw=get_fasta_writer($output,50);
while(my $f=$fr->())
{
    next unless length($f->{seq})>=$minleng;
    $fw->($f);
}




