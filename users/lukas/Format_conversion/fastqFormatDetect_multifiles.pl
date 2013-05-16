#!/usr/bin/perl
# Author: Martin Dahlo
# lukas endler: changed to take expression for list of files
# chomped each line to get rid of trailing newline
# added unidentified list of fq (in overlapping region of qualitycodes)
# Usage:  perl scriptname.pl <infiles>
# ex.
# perl scriptname.pl reads*.fq

use warnings;
use strict;


=pod

Used to detect the format of fastq files. In its current state,
it can only differentiate between sanger and solexa/illumina.
If need arises, checking for different versions of illumina formats
could easily be implemented. ( Please upload an update if you implement this )

Can easily be copy/pasted into any other script and altered to do other
things than die when it has determined the format.

Pseudo code

* Get a list of fastq files
* Look at each quality ASCII char and convert it to a number
* Depending on if that number is above or below certain thresholds,
  determine the format.


=cut


# get variables
my $usage = "Usage:  perl scriptname.pl <infiles >\n eg.: perl .\/fastqFormatDetect_multifiles.pl \'*.fq\'\n";
my $fqs = shift or die $usage;
($fqs =~ /^-+h/) && die $usage;
# get list of files
my @filelist = glob($fqs);
my @solexa = '';
my @sanger = '';
my @not_def = '';

foreach my $fq (@filelist){
  open FQ, "<", $fq or die $!;
  # initiate
  my @line;
  my $l;
  my $number;
  my $type; 
  # go thorugh the file
  while(<FQ>){
    $type = undef;
    # if it is the line before the quality line
    if($_ =~ /^\+/){      
      $l = <FQ>; # get the quality line
      chomp($l);
      @line = split(//,$l); # divide in chars
      for(my $i = 0; $i <= $#line; $i++){ # for each char
	$number = ord($line[$i]); # get the number represented by the ascii char
	
	# check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
	if($number > 74){ # if solexa/illumina and not illumina 1.8
	  # "This file is solexa/illumina format\n"; # print result to terminal and die
	  push(@solexa,$fq); $type = 'sol/il'; last;
	}elsif($number < 59){ # if sanger
	  #die "This file is sanger format\n"; # print result to terminal and die
	  push(@sanger,$fq);  $type = 'sanger';last;
	}
      }
     ($type) && last; 
    }
  }
  print $fq,"\n";
  close FQ;
  ($type) || push(@not_def,$fq);
}
print "solexa/illumina: ",join("\t",@solexa),"\n";
print "sanger: ",join("\t",@sanger),"\n";
print "no identified: ",join("\t",@not_def),"\n";
