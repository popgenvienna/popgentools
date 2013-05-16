#!/usr/bin/perl
# author: Lukas Endler Tue Dec 18 14:00:07 CET 2012
# takes an inputstream in syncfile format and calculates average coverage and stdv chromosome wise and overall for each population
# useage: perl get_coverage_from_sync.pl < infile.sync > coverage.txt
use warnings;
use strict;
#use Data::Dumper;
#use Getopt::Long;
#use Pod::Usage;


my @chrom_levs = qw( Sum SSum Zeros Length Mean Stdv );
my %chroms=();
my $count = 0;
my $pops = 0;
my $spaces = qr/\s+/o;
my $lastchrom = "";

if (-t STDIN and not @ARGV) {
  die "useage: perl get_coverage_from_sync.pl < infile.sync > coverage.txt\ntakes an inputstream in syncfile format and calculates average coverage and stdv chromosome wise and overall for each population\n";  
} 
print STDERR "reading sync file:\n"; 
while (<>){
  $count += 1;
  #my $line = $_;
  chomp($_);
  #my @fields = split(/\s+/, $line);
  my @fields = split($spaces, $_);
  my $c_name = $fields[0];
  if ($c_name ne $lastchrom) {
    $lastchrom = $c_name;
    print STDERR "Chrom $c_name : 1 \n";
  }
  unless (exists $chroms{$c_name}) {
    $chroms{$c_name} = ();
    for my $i (qw( Sum SSum Zeros)) {$chroms{$c_name}{$i} = [ ] };
    for my $i (qw( Length Mean Stdv)) {$chroms{$c_name}{$i} = 0 };
    #   @chroms{$c_name}{ qw( Sum SSum Zeros) } = ([ ], [ ], [ ]); 
    #   @chroms{$c_name}{ qw( Length Mean Stdv) } = (0,0,0);
  }
  $chroms{$c_name}{'Length'} += 1;
  unless ($chroms{$c_name}{'Length'} % 1e5) {print STDERR "Chrom $c_name : ",$chroms{$c_name}{'Length'},"\n"}; 
  $pops = $#fields -3;
  for my $i (3..$#fields) {
    my $val = 0; 
    for my $j (split(":",$fields[$i])) { $val += $j};
    unless (exists $chroms{$c_name}{'Sum'}[$i-3]) {
      for (qw( Sum SSum Zeros)) {${$chroms{$c_name}{$_}}[$i-3] = 0}; 
      #push(@{$chroms{$c_name}{'Sum'}}, 0);
    }
    $chroms{$c_name}{'Sum'}[$i-3] += $val;
    $chroms{$c_name}{'SSum'}[$i-3] += $val*$val;
  #  if ($val == 0) {$chroms{$c_name}{'Zeros'}->[$i-3] += 1 };
  }
 # ($count == 1000) && last;
} 

print STDERR "finished reading, calculating coverage\n"; 

for my $key (keys(%chroms)){
  $chroms{$key}{'Mean'} = mean($chroms{$key}{'Sum'},$chroms{$key}{'Length'});
  $chroms{$key}{'Stdv'} = stdv($chroms{$key}{'Sum'},$chroms{$key}{'SSum'},$chroms{$key}{'Length'});
}


print "chrom\t",join("\t",map(sprintf("Pop%i",$_+1), (0..$pops))),"\tLength\n"; 

for my $key ( keys(%chroms) ){ 
  print $key; 
  for my $j (0..$#{$chroms{$key}{'Mean'}}){
	printf("\t%.1d %.1d",$chroms{$key}{'Mean'}[$j],$chroms{$key}{'Stdv'}[$j])    
    }   
  printf("\t%i\n",$chroms{$key}{'Length'});
}
## Calculate totals
## tot. length
my $overall_length = 0;
$overall_length += $chroms{$_}{'Length'} for ( keys(%chroms) );
## total mean
my @overall_mean = map {my $i=0; for my $key (keys(%chroms)) {$i += $chroms{$key}{'Mean'}[$_]*$chroms{$key}{'Length'}/$overall_length }; $i} (0..$pops);
## total stdv
my @overall_stdv =  map {my $i=0; for my $key (keys(%chroms)) {$i += $chroms{$key}{'SSum'}[$_]}; $i = sqrt($i/$overall_length - $overall_mean[$_]**2)} (0..$pops);
## Print Totals
print "Total" ;
for my $j (0..$pops){
  printf("\t%.1d %.1d",$overall_mean[$j],$overall_stdv[$j])    
}
printf("\t%i\n",$overall_length);
my @overall_str = [];
for my $j (0..$#overall_mean) {$overall_str[$j]=sprintf("%.1d",$overall_mean[$j])};
print "Total\t",join(",",@overall_str),"\n";

#######################
### end of main
#######################
sub mean {
  (my $sums , my $length) = @_;
  my @means = map { $_ / $length }  @{$sums};
  return \@means ; 
}

sub stdv {
  (my $sums, my $ssums, my $length ) = @_;
  my @stdvs = map {sqrt ( (${$ssums}[$_] - (${$sums}[$_]**2)/$length)/$length  )} (0..$#{$sums});
  return \@stdvs; 
}

