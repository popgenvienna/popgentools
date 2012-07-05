#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;


my $SYNC_FILE_NAME;
my @COLUMNS=();
my @SET_MAJOR=();
my $HELP;

GetOptions(
	"input=s"=>\$SYNC_FILE_NAME,
	"populations=s"=> \@COLUMNS,
	"set-major-according-to=s"=> \@SET_MAJOR,
	'help'=> \$HELP,
) or die "GetOptions";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

@COLUMNS = (split(/,/, join(',', @COLUMNS)));
@SET_MAJOR = (split(/,/, join(',', @SET_MAJOR)));

if (scalar @SET_MAJOR == 0){die "Option --set-major-according-to not specified."}
if (scalar @COLUMNS == 0){die "Option --populations not specified."}



sub max{
	my $m=$_[0];
	for (my $i=0; $i< scalar(@_); $i++){
		next unless ($_[$i]>$m);
		$m = $_[$i];		
	} 
	return $m;
}


my %major=();
foreach my $e (@SET_MAJOR){
	$major{$e}=1;	
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
			
			if (defined($major{$COLUMNS[$i]})){
				$sum_a+=$a;
				$sum_t+=$t;
				$sum_c+=$c;
				$sum_g+=$g;
			}
		}
		
		my $m = max($sum_a, $sum_t, $sum_c, $sum_g);
			
		my $base="";
		if ($m==$sum_a){	$base="a";}
		elsif($m==$sum_t){	$base="t";}
		elsif($m==$sum_c){	$base="c";}
		else{			$base="g";}
		
		for (my $i=0; $i<= (scalar(@COLUMNS) -1); $i++ ){
			my $locMax=0; 
			if ($base eq "a"){	$locMax= $populations{$COLUMNS[$i]-1}{'a'};}
			elsif($base eq "t"){	$locMax= $populations{$COLUMNS[$i]-1}{'t'};}
			elsif($base eq "c"){	$locMax= $populations{$COLUMNS[$i]-1}{'c'};}
			else{					$locMax= $populations{$COLUMNS[$i]-1}{'g'};}
			
			my $f; 
			
			my $sum_per_SNP = ($populations{$COLUMNS[$i]-1}{'a'} + $populations{$COLUMNS[$i]-1}{'t'} + $populations{$COLUMNS[$i]-1}{'c'} + $populations{$COLUMNS[$i]-1}{'g'});
			
			if ($sum_per_SNP!=0){
			if ($m!=0){ $f = $locMax/$sum_per_SNP; }else{$f=0;}
			}else{
				$f=0;	
			}
			$line = $line."\t".$f;
		}
		
		$line = $line."\n";
		print $line;	
	}
}
close fileHandle;

1;

__END__

=pod

=head1 NAME

add-pop-major-frequencies.pl

=head1 SYNOPSIS

perl add-pop-major-frequencies.pl --input inputFileName --populations column1,column2,column3... --set-major-according-to column1,column2,column3...

=head1 DESCRIPTION

The script takes an input file (for example synchronized pileup file or output of CMH test) and calculates a major allele frequency 
for each of the specified populations (by option '--populations').

For each position of the input file the script prints whole input line and also additional major allele frequencies 
of the specified populations. The frequencies are ordered in the same way as population columns specified by option '--populations'. 
A major allele is the allele that is major over all populations specified by option '--set-major-according-to'.

=head1 OPTIONS

=over 4

=item --input

A name of an input file. The file should be a tab delimited table that contains several populational columns in the same format as in synchronized pileup file. For example, a file in synchronized pileup format or output of CMH file can be used as an input.
 
=item --populations

A comma separated list containing columns of populations for which the frequencies shouls be calculated, numbering of columns is 1 based (the first column has number 1). Usually, the first population is in column 4. At least one column should be specified by this option.

=item --set-major-according-to

A comma separated list containing columns of populations that will be used to set major allele, numbering of columns is 1 based (the first column has number 1). Usually, the first population is in column 4. At least one column should be specified by this option.

=back

=head1 EXAMPLE

=over 4

=item INPUT file example.input.sync

	3L	7	T	0:112:13:0:0:0	0:97:5:0:0:0	0:108:8:0:0:0
	3L	9	C	0:20:100:0:0:0	0:16:102:0:0:0	0:120:15:0:0:0:0

=item USAGE

perl add-pop-major-frequencies.pl --input example.input.sync --populations 4,5,6 --set-major-according-to 4,5 > example.output.freq

=item OUTPUT file example.output.freq

	3L	7	T	0:112:13:0:0:0	0:97:5:0:0:0	0:108:8:0:0:0	0.896	0.950980392156863	0.931034482758621
	3L	9	C	0:20:100:0:0:0	0:16:102:0:0:0	0:120:15:0:0:0:0	0.833333333333333	0.864406779661017	0.111111111111111
	
=back

=cut
