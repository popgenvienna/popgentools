#!/usr/bin/perl -w
use strict;
use warnings;

my $pa=[];

    #@HWUSI-EAS300R:6:1:0:1361#0/1
    #NCGACGGCACCCACGAGCTGTGGGGATGGGGCGACAAGGATCTCACCGAGATCGGAAGAGCGGTTCAGCAGGAATG


my $flag=0;
my $counter=0;
while(<>)
{
    if($flag==0 && $_=~/^[+]/)
    {
         $flag=1;
         next;
    }
    next unless $flag;
    
    chomp;
    my @a=split//,$_;
    
    for(my $i=0;$i<@a;$i++)
    {
       my $val= ord($a[$i]) - 64;
       $pa->[$i]+=$val;
    }
    
    
    $counter++;
    $flag=0;
}


for(my $i=0;$i<@$pa; $i++)
{
    my $val=$pa->[$i];
    $val/=$counter;
    $val=sprintf("%.2f",$val);
    my $printpos=$i+1;
    print "$printpos\t$val\n";
    
}


