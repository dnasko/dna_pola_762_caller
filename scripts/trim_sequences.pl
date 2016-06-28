#!/usr/bin/perl -w

# MANUAL FOR trim_sequences.pl

=pod

=head1 NAME

trim_sequences.pl -- Given the out file from 762_caller, will trim PolAs

=head1 SYNOPSIS

 trim_sequences.pl --in=/Path/to/infile.fasta --trim=/Path/to/trim_file.txt --out=/Path/to/output_trimmed.fasta --only_complete
                     [--help] [--manual]

=head1 DESCRIPTION

 Given the output file from 762_caller and the FASTA of DNA polymerase A's
 this script will trim each sequence so that only the 547-926 region exists
 for each sequence.

 With the --only_complete flag this parser will only print and trim sequences
 that are complete (i.e. span, end-to-end the 547-926 region).
 
=head1 OPTIONS

=over 3

=item B<-i, --in>=FILENAME

Input file in peptide FASTA format. (Required) 

=item B<-t, --trim>=FILENAME

Input file from the output of 762_caller. This is a tab-delimmited file that deatils the 547 - 926 trimming sites. (Required)

=item B<-o, --out>=FILENAME

Output file in FASTA format. (Required) 

=item B<-c, --only_complete>

Using this flag the parser will only output sequences that spanned, end-to-end, the 547-926 region. Those that do not are thrown out. (Optional)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries.

Bio::SeqIO

=head1 AUTHOR

Written by Daniel Nasko, 
Center for Bioinformatics and Computational Biology, University of Delaware.

=head1 REPORTING BUGS

Report bugs to dnasko@udel.edu

=head1 COPYRIGHT

Copyright 2016 Daniel Nasko.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's 
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut


use strict;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use FindBin;
use Bio::SeqIO;

#ARGUMENTS WITH NO DEFAULT
my($infile,$trim,$outfile,$only_complete,$help,$manual);
my $version = "1.0";
GetOptions (
                           "i|in=s"          =>\$infile,
                           "t|trim=s"        =>\$trim,
                           "o|out=s"         =>\$outfile,
                           "c|only_complete" =>\$only_complete,
                           "h|help"          =>\$help,
                           "m|manual"        =>\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --in not found.\n\n", -exitval => 2, -verbose => 1)  if (! $infile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --trim not found.\n\n", -exitval => 2, -verbose => 1)  if (! $trim);

my %Coords;
my $l=0;
open(IN,"<$trim") || die "\n Cannot open the trimming file: $trim\n";
while(<IN>) {
    chomp;
    if ($l>0) { ## If not on the first line
	my @a = split(/\t/, $_);
	$Coords{$a[0]}{"start"} = $a[3];
	$Coords{$a[0]}{"stop"}  = $a[4];
	$Coords{$a[0]}{"complete"} = $a[5];
    }
    $l++;
}
close(IN);

my $seq_in  = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => $infile );

open(OUT,">$outfile") || die "\n Cannot write to the output file: $outfile\n";
while( my $seq = $seq_in->next_seq() ) {
    my $header = $seq->id;
    my $sequence = $seq->seq;
    if (exists $Coords{$header}) {
	unless ($only_complete && $Coords{$header}{"complete"} eq "no") { ## this is the magic of the --only_complete flag
	    my $left = $Coords{$header}{"start"};
	    my $right = $Coords{$header}{"stop"} - $Coords{$header}{"start"};
	    my $trim = substr $sequence, $left, $right;
	    print OUT ">" . $header . "\n" . $trim . "\n";
	}
    }
    else { die "\n Cannot find $header in the trim file! Make sure 762_caller was run!\n\n"; }    
}
close(OUT);

exit 0;
