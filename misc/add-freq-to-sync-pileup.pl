#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;


my $SYNC_FILE_NAME;
my @COLUMNS=();
my $HELP = 0;

GetOptions(
	"sync-file=s"=>\$SYNC_FILE_NAME,
	'columns=i{,}'=> \@COLUMNS,
	'help'=> \$HELP,
) or die "GetOptions";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DETAILS"}) if $HELP;

sub max{
	my $m=$_[0];
	for (my $i=0; $i< scalar(@_); $i++){
		next unless ($_[$i]>$m);
		$m = $_[$i];		
	} 
	return $m;
}

open fileHandle, "<", $SYNC_FILE_NAME or die "Could not open file $SYNC_FILE_NAME";
while (my $line = <fileHandle>){
	if ($line =~ m/^#/){
		print $line;
	}else{
		chomp($line);
		my @parts = split "\t", $line;		

		my %populations=();
		
		my $sum_a=0; my $sum_t=0; my $sum_c=0; my $sum_g=0;
		
		for (my $i=0; $i<= (scalar(@COLUMNS) -1); $i++ ){
			
			my ($a,$t,$c,$g) = split ":", $parts[$COLUMNS[$i]-1];			
			
			$populations{$COLUMNS[$i]-1}{'a'} = $a;
			$populations{$COLUMNS[$i]-1}{'t'} = $t;
			$populations{$COLUMNS[$i]-1}{'c'} = $c;
			$populations{$COLUMNS[$i]-1}{'g'} = $g;			
			
			$sum_a+=$a;
			$sum_t+=$t;
			$sum_c+=$c;
			$sum_g+=$g;
		}
		
#		print Dumper(\%populations);
		
		my $m = max($sum_a, $sum_t, $sum_c, $sum_g);
		
#		print Dumper($m);
			
		my $base="";
		if ($m==$sum_a){	$base="a";}
		elsif($m==$sum_t){	$base="t";}
		elsif($m==$sum_c){	$base="c";}
		else{			$base="g";}
		
#		print Dumper($base);
		
		for (my $i=0; $i<= (scalar(@COLUMNS) -1); $i++ ){
			my $locMax=0; 
			if ($base eq "a"){	$locMax= $populations{$COLUMNS[$i]-1}{'a'};}
			elsif($base eq "t"){	$locMax= $populations{$COLUMNS[$i]-1}{'t'};}
			elsif($base eq "c"){	$locMax= $populations{$COLUMNS[$i]-1}{'c'};}
			else{					$locMax= $populations{$COLUMNS[$i]-1}{'g'};}
			
			my $f; 
			if ($m!=0){ $f = $locMax/ ($populations{$COLUMNS[$i]-1}{'a'} + $populations{$COLUMNS[$i]-1}{'t'} + $populations{$COLUMNS[$i]-1}{'c'} + $populations{$COLUMNS[$i]-1}{'g'}); }else{$f=0;}
			$line = $line."\t".$f;
		}
		
		$line = $line."\n";
				
		#my $col = $COLUMNS[0];
		#my ($a,$t,$c,$g) = split ":", $parts[$col-1];
		#my $m = max($a,$t,$c,$g);
		#my $f;
		#if ($m!=0){ $f = $m/($a+$t+$c+$g);}else{ $f = 0;}
		#$line = $line."\t".$f;
			
		#my $base="";
		#if ($m == $a){
		#	$base="a";	
		#}elsif($m==$t){
		#	$base="t";		
		#}elsif($m==$c){
		#	$base="c";		
		#}else{
		#	$base="g";		
		#}
		
		
		#for (my $i=1; $i<= (scalar(@COLUMNS) -1); $i++ ){

		#	my ($a,$t,$c,$g) = split ":", $parts[$COLUMNS[$i]-1];			

		#	if ($base eq "a"){
		#		$m=$a;
		#	}elsif($base eq "t"){
		#		$m=$t;
		#	}elsif($base eq "c"){
		#		$m=$c;
		#	}else{$base=$g}
		#	
		#	my $f;
		#	if ($m!=0){ $f = $m/($a+$t+$c+$g);}else{ $f = 0;}
		#	$line = $line."\t".$f;	
		#}
		print $line;	
	}
}
close fileHandle;

1;

__END__

=pod

=head1 NAME

add-freq-to-sync-pileup.pl

=head1 SYNOPSIS

perl add-freq-to-sync-pileup.pl --sync-file synchronizedPileupFileName --columns column1 column2 column3 ...

=head1 DETAILS

The script adds new columns to synchronized pileup file. The columns contain for each population (in the same order as populations in the sync-pileup) frequencies of the major allele. 

For population with counts a:t:c:g:n:i, where a>t,c,g, the frequncy a/(a+t+c+g) will be added.

Columns are numbered from 1. At least one have to be specified.

=cut
