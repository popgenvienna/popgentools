#!/usr/bin/perl -w
use strict;



my $totest=100000;

my $coverage=100;
my $snpcount=2;

my $rawsample=$snpcount x "1" . ($coverage-$snpcount) x "0";
my $sample=split //,$rawsample

my $totest=$

my $freqs=[];
for my $i (1..$totest)
{
    my $states=[];
    for my $k (1..$coverage)
    {
        my $index=int(rand($sample));
        push @$states,$sample->[$index];
    }
    my $index=int(rand()*$sample);
}

my $values=[];

