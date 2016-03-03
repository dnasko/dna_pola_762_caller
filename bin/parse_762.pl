#!/usr/bin/perl
use strict;

my $infile = $ARGV[0];
my $line = 0;
my $array_pos = "";

open(IN,"<$infile");
while(<IN>) {
    chomp;
    my $pos = 0;
    if ($line == 1) {
	my @A = split(//, $_);
	for(my $i=0;$i<scalar(@A); $i++) {
	    unless ($A[$i] eq "-") {
		$pos++;
		if ($pos == 762) {
		    $array_pos = $i;
		}
	    }
	}
    }
    elsif ($line == 6) { 
	my $header = $_;
	$header =~ s/^>//;
	    print "$header\t";
    }
    elsif ($line == 7) {
	my @A = split(//, $_);
	my $local_pos = 0;
	for(my $i=0;$i<=$array_pos; $i++) {
	    unless ($A[$i] eq "-") {
		$local_pos++;
	    }
	}
	print $A[$array_pos] . "\t" . $local_pos . "\n";
    }
    $line++;
}
close(IN);
