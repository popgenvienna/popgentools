#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use List::Util qw[min max];
use TEHierarchy;
use ParseSam;
use TEInsertUtility;
use TESamReader;


my $input=shift;
my $tehierfile=shift;

die "usage: $0 input tehierfile" unless $input;


my $teh_resolver=get_te_hierarchy_resolver($tehierfile,"family");
my $samparser=get_te_samparser($teh_resolver);
my $sr=TESamReader->new($input,$samparser,1000000000000000);


my($countabs,$countpres)=(0,0);

while(my $s=$sr->next())
{
    if($s->{mode} eq "abs")
    {
        $countabs++;

    }
    elsif(  $s->{mode} eq "pre")
    {
        $countpres++;

    }
    else
    {
        die "Mode not supported";
    }
}


my $sum=$countabs+$countpres;
print "potential useful pairs: $sum\n";
print "potential absence pairs: $countabs\n";
print "potential presence pairs: $countpres\n";
