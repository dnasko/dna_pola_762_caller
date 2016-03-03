#!/usr/bin/perl -w

# MANUAL FOR find_762.pl

=pod

=head1 NAME

find_762.pl -- Iteratively runs MAFFT alignments against a set of reference(s).

=head1 SYNOPSIS

 find_762.pl --in=/Path/to/infile.fasta --ref=/Path/to/references.fasta --out=/Path/to/output.txt
                     [--help] [--manual]

=head1 DESCRIPTION

 This script runs a MAFFT alignment for each sequence in your input.fasta
 file against all of the refrecne sequecnes in references.fasta.
 
=head1 OPTIONS

=over 3

=item B<-i, --in>=FILENAME

Input file in peptide FASTA format. (Required) 

=item B<-r, --ref>=FILENAME

Input set of peptide reference sequecnes. (Required)

=item B<-o, --out>=FILENAME

Output file in text format. (Required) 

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries.
Bio Perl

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
my($infile,$ref,$outfile,$help,$manual);

GetOptions (
    "i|in=s"=>\$infile,
    "r|ref=s"=>\$ref,
    "o|out=s"=>\$outfile,
    "h|help"=>\$help,
    "m|manual"=>\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --in not found.\n\n", -exitval => 2, -verbose => 1)  if (! $infile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --ref not found.\n\n", -exitval => 2, -verbose => 1)  if (! $ref);

## Create the temporary working directory.
my @chars = ("A".."Z", "a".."z");
my $rand_string;
$rand_string .= $chars[rand @chars] for 1..8;
my $working_dir = "mega_mafft_tmp_" . $rand_string;
print `mkdir $working_dir`;

## Begin
print " Beginning Mega MAFFT: ";
print `date`;

my $seq_in  = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => $infile );

if (-e $outfile) { print `rm $outfile`; } ## Clear out the outfile in case it's still there.
while( my $seq = $seq_in->next_seq() ) {
    # print `cat $FindBin::Bin/../refs/references.fasta > $outdir/tmp.fa`;
    print `cat $ref > $working_dir/tmp.fa`;
    open(OUT,">>$working_dir/tmp.fa") || die "Can't open the temporary FSA: $working_dir/tmp.fa\n\n";
    print OUT ">" . $seq->id . "\n" . $seq->seq . "\n";
    close(OUT);
    print `mafft --retree 2 --inputorder $working_dir/tmp.fa > $working_dir/aln.fa 2> $working_dir/std.err`;
    print `$FindBin::Bin/flaten_fasta.pl -i $working_dir/aln.fa > $working_dir/aln.flt.fa`;
    # print `$FindBin::Bin/parse_762.pl $working_dir/aln.flt.fa >> $working_dir/results.txt`;
    print `$FindBin::Bin/parse_762.pl $working_dir/aln.flt.fa >> $outfile`;
}

print "Done Mega MAFFT: ";
print `date`;

print `rm -rf $working_dir`;

exit 0;
