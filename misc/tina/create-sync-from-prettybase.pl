#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $HELP;

my $prettybaseFileName;
my $chromosome;
my $type=2;

GetOptions(
	"chromosome=s"=>\$chromosome,
	"prettybase=s"=>\$prettybaseFileName,
	"type=i"=>\$type, #(single=1/double=2)
	 'help'=> \$HELP,
);

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

if (!defined($prettybaseFileName) or !defined($chromosome) or ($type <1) or ($type>2)){die "Invalid arguments, use 'perl $0 --help'."}

sub max{
	my $max = $_[0];
	
	foreach my $value (@_){
		if ($value > $max){$max=$value;}
	}
	return $max;
}


my $ptrPositions={};
my $last=0;

srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);

open fileHandlePrettybase, "<", $prettybaseFileName or die "Could not open file $prettybaseFileName";

while (my $line = <fileHandlePrettybase>){
	next if $line =~ m/^#/;
	chomp($line);
	my ($position, $flyName, $base1, $base2) = split "\t", $line;	
	$position = int($position);

	my $r=rand(2);
	
	if (defined($ptrPositions->{$position})){
		if ($type == 1){
		
			if ($r==0){$ptrPositions->{$position}{$base1}+=1;}
			else{$ptrPositions->{$position}{$base2}+=1;}
		}else{
			$ptrPositions->{$position}{$base1}+=1;
			$ptrPositions->{$position}{$base2}+=1;
		}
	}else{
		$ptrPositions->{$position}{A}=0;
		$ptrPositions->{$position}{T}=0;	
		$ptrPositions->{$position}{C}=0;
		$ptrPositions->{$position}{G}=0;
		$ptrPositions->{$position}{N}=0;
		
		if ($type == 1){		
			if ($r==0){$ptrPositions->{$position}{$base1}+=1;}
			else{$ptrPositions->{$position}{$base2}+=1;}
		}elsif($type == 2){
			$ptrPositions->{$position}{$base1}+=1;
			$ptrPositions->{$position}{$base2}+=1;
		}
	}
	if ($last<$position){$last=$position;}
}
close fileHandlePrettybase;

#print Dumper($ptrPositions);
#print Dumper($ptrPositions->{184});

print "#chromosome\tposition\trefBase\tpopulation\n";
for(my $position = 0; $position<=$last; $position++){

	next if !(defined($ptrPositions->{$position})); 

	if (!(defined($ptrPositions->{$position}{A}))){
		$ptrPositions->{$position}{A} = 0;
	}
	if (!(defined($ptrPositions->{$position}{T}))){
		$ptrPositions->{$position}{T} = 0;
	}
	if (!(defined($ptrPositions->{$position}{C}))){
		$ptrPositions->{$position}{C} = 0;
	}
	if (!(defined($ptrPositions->{$position}{G}))){
		$ptrPositions->{$position}{G} = 0;
	}
	if (!(defined($ptrPositions->{$position}{N}))){
		$ptrPositions->{$position}{N} = 0;
	}
	if (!(defined($ptrPositions->{$position}{X}))){
		$ptrPositions->{$position}{X} = 0;
	}
	
	my $A = $ptrPositions->{$position}{A};
	my $T =	$ptrPositions->{$position}{T};
	my $C = $ptrPositions->{$position}{C};
	my $G = $ptrPositions->{$position}{G};
	my $N = $ptrPositions->{$position}{N} + $ptrPositions->{$position}{X};
	my $X = 0;
		
	foreach my $key (keys %{$ptrPositions->{$position}}){
		next if (($key eq "A")||($key eq "T")||($key eq "C")||($key eq "G")||($key eq "N")||($key eq "X"));
		$X+= $ptrPositions->{$position}{$key};			
	}
	
		
	my $max = max($A, $T, $C, $G, $N, $X);
	
	my $refBase;

	if ($max == $N){$refBase = "N"};
	if ($max == $X){$refBase = "N"};  	
	if ($max == $A){$refBase = "A"};
	if ($max == $T){$refBase = "T"};
	if ($max == $C){$refBase = "C"};
	if ($max == $G){$refBase = "G"};
	
	print $chromosome."\t".
		  $position."\t".
		  $refBase."\t".
		  $A.":".$T.":".$C.":".$G.":".$N.":".$X."\n";	
}

=pod

=head1 NAME

create-sync-from-prettybase.pl

=head1 SYNOPSIS

perl create-sync-from-prettybase.pl --prettybase fileName --type 1(2) --chromosome 3L

=head1 DESCRIPTION
 
The script takes prettybase data and creates synchronized pileup out of it. For details about synchronized pileup file format see popoolation script synchronize-pileup.pl. For details about prettybase, see google.

=head1 OPTIONS 

=over 4

=item --prettybase

An input file in a prettybase file format.

=item --chromosome 

A name of a chromosome for which the prettybase file was created. 

=item --type 

Choose weather to use 2 alleles data for each position or only 1 allele (randomly choosed for each position). 
Possible values fot the parameter are 1 and 2, default value is 2.

=item --help, -h

Prints the help page.

=back

=head1 EXAMPLE

=over 4

=item INPUT example.prettybase

  1	id1	A	T
  1	id2	T	T
  1	id3	A	A
  2	id1	G	C
  2	id2	C	C
  2	id3	G	C
  3	id1	A	A
  3	id2	A	T
  3	id3	T	A

=item Used command

perl create-sync-from-prettybase.pl --prettybase example.prettybase --type 2 --chromosome 3L

=item OUTPUT

  #chromosome	position	refBase	population
  3L	1	T	3:3:0:0:0:0
  3L	2	C	0:0:4:2:0:0
  3L	3	A	4:2:0:0:0:0
