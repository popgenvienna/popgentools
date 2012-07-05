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
my $tosample=shift;
my $insample=shift;

die "usage: $0 input tehierfile tosample insample" unless $input;


my $teh_resolver=get_te_hierarchy_resolver($tehierfile,"family");
my $samparser=get_te_samparser($teh_resolver);
my $samformater=get_basic_samformater
my $sr=TESamReader->new($input,$samparser,1000000000000000);



my $ratio;

# parse the file once and print already the putative absence fragments;
# store the presence fragments in a hash, need to be parsed at the second go
my $readhash={};
while(my $s=$sr->next())
{
    $ratio=$tosample/$insample;
    my $rand=rand();
    
    if($rand<$ratio)
    {
        if($s->{mode} eq "pre")
        {
            my $readid=$s->{r1}{readid};
            $readhash->{$readid}=1;
        }elsif($s->{mode} eq "abs")
        {
            print $samformater->($s->{r1})."\n";
            print $samformater->($s->{r2})."\n";
        }
     # print
     $tosample--;
    }
    $insample--;
}
$sr->close();

open my $ifh, "<", $input or die "Could not open input file";
while(my $line=<$ifh>)
{
    chomp $line;
    my($readid,undef)=split /\t/,$line;
    if(exists($readhash->{$readid}))
    {
        print $line."\n";
    }
}






