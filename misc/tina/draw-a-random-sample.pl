#!/usr/bin/perl


use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


my $file;
my $count;
my $HELP;

GetOptions(
	"input=s"=>\$file,
	"size=i"=>\$count,
	"help"=>\$HELP
)or die "Invalid arguments, use 'perl $0 --help'.";

pod2usage({-verbose=>99, -sections=>"NAME|SYNOPSIS|DESCRIPTION|OPTIONS|EXAMPLE"}) if $HELP;

if (!defined($file)){die "Option --input should be specified."}

my $text = `wc -l $file`;

my $total= (split(" ", $text))[0];

srand (time ^ $$ ^ unpack "%L*", `ps axww | gzip -f`);

open fileHandle, "<", $file; 

my $line = <fileHandle>;
my $lastAll = $total;
my $lastToDraw = $count;
my $index=0;


while (($lastToDraw > 0) and ($index <= $total)){

	my $r = rand($lastAll);
	if ($r < $lastToDraw){
		print $line;	
		$lastToDraw-=1;	
	}
	$lastAll-=1;
	$index+=1;
	
	$line = <fileHandle>;
}



close fileHandle;

=pod

=head1 NAME 

draw-a-random-sample.pl

=head1 SYNOPSIS

perl draw-a-random-sample.pl --input fileName --size N

=head1 DESCRIPTION 

The script randomly draw a specified number of lines form an input file and prints them out. Each line is drawn with equal probability and the order of drawn lines does not change during sampling process.

=head1 OPTIONS

=over 4

=item --input

An input file name. The file can be of any type. 

=item --size

A number of lines that should be drawn. 

=item --help, -h

Prints the help page.

=back



=head1 EXAMPLE OF USAGE

=over 4

=item INPUT file example.list

	1
	2
	3
	4
	5
	6
	7

=item Used command

perl draw-a-random-sample.pl --input example.list --size 3 > output.list

=item OUTPUT file output.list (this is only example, your output is probably different)

	1
	5
	6

=back

=cut
