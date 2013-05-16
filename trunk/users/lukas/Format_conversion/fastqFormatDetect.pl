#!/usr/bin/perl
# Author: Martin Dahlo
# changed by Lukas Endler: chomp line so newline does not gie false sanger encoding
# changed by Lukas Endler: changed upper threshold to 74, so illumina 1.8+ format comes as sanger
# changed by Lukas Endler: tell in case not determinable
# Usage:  perl scriptname.pl <infile>
# ex.
# perl scriptname.pl reads.fq

use warnings;
use strict;


=pod

Used to detect the format of a fastq file. In its current state,
it can only differentiate between sanger and solexa/illumina.
If need arises, checking for different versions of illumina formats
could easily be implemented. ( Please upload an update if you implement this )

Can easily be copy/pasted into any other script and altered to do other
things than die when it has determined the format.

Pseudo code

* Open the fastq file
* Look at each quality ASCII char and convert it to a number
* Depending on if that number is above or below certain thresholds,
  determine the format.


=cut


# get variables
my $usage = "Usage:  perl scriptname.pl <infile >\n";
my $fq = shift or die $usage;

# open the files
open FQ, "<", $fq or die $!;


# initiate
my @line;
my $l;
my $number;


# go thorugh the file
while(<FQ>){
  
  # if it is the line before the quality line
  if($_ =~ /^\+/){
    $l = <FQ>; # get the quality line
    chomp($l); #get rid of trailing newline
    @line = split(//,$l); # divide in chars
    for(my $i = 0; $i <= $#line; $i++){ # for each char
      $number = ord($line[$i]); # get the number represented by the ascii char
      # check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
      # actually should be lower than 73, but new Illumina 1.8 which is similar to Sanger goes to 74
      if($number > 74){ # if solexa/illumina not illumina 1.8
        die "This file is solexa/illumina format\n"; # print result to terminal and die
      }elsif($number < 59){ # if sanger
        die "This file is sanger format\n"; # print result to terminal and die
      }
    }
  }
}
die "Could not determine the format of quality encoding\n"
