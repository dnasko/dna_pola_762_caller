#!/usr/bin/perl

#MANUAL FOR flaten_fasta.pl

=pod

=head1 NAME

 flaten_fasta.pl -- flatens a FASTA file

=head1 SYNOPSIS

 flaten_fasta.pl -in /Path/to/infile.fasta -out /Path/to/outfile.fasta
                         [-help]  [-manual]

=head1 DESCRIPTION

 Flatens a FASTA file, by this I mean will remove all return characters embedded
 in the sequence, such that when this program is done running the number of
 sequences will be equal to half the number of lines.

=head1 OPTIONS

=over 3

=item B<-i, --in>=FILENAME

Input file in FASTA format. (Required)

=item B<-o, --out>=FILENAME

Output file in FASTA format. (Optional)
 default prints to STDOUT.

=item B<-h, --help>

Displays the usage message. (Optional)

=item B<-m, --manual>

Displays fill manual. (Optional)

=back

=head1 DEPENDENCIES

Requires the following Perl libraries.

-none-

=head1 AUTHOR

Written by Daniel Nasko,
Center for Bioinformatics and Computational Biology, University of Delaware.

=head1 REPORTING BUGS

Report bugs to dnasko@udel.edu

=head1 COPYRIGHT

Copyright 2013 Daniel Nasko.
License GPLv3+: GPU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRENTY, to the extent permitted by law.

Please acknowledge author and affliliation in published work arising from thie script's
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut

use strict;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

## ARGUMENTS WITH DEFAULTS
my $outfile;

## ARGUMENTS WITH NO DEFAULT
my($infile,$help,$manual);

GetOptions (
                       "i|in=s"       =>     \$infile,
                       "o|out=s"      =>     \$outfile,
                       "h|help=s"     =>     \$help,
                       "m|manual=s"   =>     \$manual);

## VALIDATE ARGS
pod2usage(-verbose => 2) if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($help);
pod2usage( -msg  => "\n\n Fatal! Required argument -in not found.\n\n", -exitval => 2, -verbose => 1) if (! $infile );

if ($infile =~ m/\.gz$/) { ## if a gzip compressed infile
    open(IN,"gunzip -c $infile |") || die "\n\n Cannot open the input file: $infile\n\n";
}
else { ## If not gzip comgressed 
    open(IN,"<$infile") || die "\n\n Cannot open the input file: $infile\n\n";
}

my $line_count = 0;

if ($outfile eq "") {
    while(<IN>) {
	chomp;
	if ($line_count == 0) {
	    unless ($_ =~ m/^>/) {
		die "\n\n Error: The infile you have provided does not appear to be in FASTA format, note the first line:\n\n$_\n\n";
	    }
	    print "$_\n";
	}
	elsif ($_ =~ m/^>/) {
	    print "\n" . $_ . "\n";
	}
	else {
	    print $_;
	}
	$line_count++;
    }
    print "\n";
}
else {
    open(OUT,">$outfile") || die "\n\n Error! Cannot open the output file: $outfile\n\n";
    while(<IN>) {
	chomp;
	if ($line_count == 0) {
	    unless ($_ =~ m/^>/) {
                die "\n\n Error: The infile you have provided does not appear to be in FASTA format, note the first line:\n\n$_\n\n";
            }
            print OUT "$_\n";
        }
        elsif ($_ =~ m/^>/) {
            print OUT "\n" . $_ . "\n";
        }
        else {
            print OUT $_;
        }
        $line_count++;
    }
    print OUT "\n";
    close(OUT);
}
close(IN);

exit 0;
