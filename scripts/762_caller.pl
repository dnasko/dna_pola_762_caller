#!/usr/bin/perl -w

# MANUAL FOR find_762.pl

=pod

=head1 NAME

762_caller.pl -- Iteratively runs MAFFT peptide alignments against a set of references to determine DNA polymerase A 762 calls.

=head1 SYNOPSIS

 762_caller.pl --in=/Path/to/infile.fasta --ref=/Path/to/references.fasta --out=/Path/to/output.txt
                     [--debug] [--help] [--manual]

=head1 DESCRIPTION

 This script runs a MAFFT alignment for each sequence in your input.fasta
 file against all of the reference sequecnes in references.fasta.
 
=head1 OPTIONS

=over 3

=item B<-i, --in>=FILENAME

Input file in *peptide* FASTA format. (Required) 

=item B<-r, --ref>=FILENAME

Input set of peptide reference sequecnes. Use the 3_references.fasta file located in the references directory. (Required)

=item B<-o, --out>=FILENAME

Output file in text format. (Required) 

=item B<-d, --debug>

Supress deleting the temporary working directory, for debugging purposes. (Optional)

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

#ARGUMENTS WITH NO DEFAULT
my($infile,$ref,$outfile,$help,$manual,$debug);
my $version = "2.0";
GetOptions (
                           "i|in=s"   =>\$infile,
                           "r|ref=s"  =>\$ref,
                           "o|out=s"  =>\$outfile,
                           "h|help"   =>\$help,
                           "m|manual" =>\$manual,
                           "d|debug"  =>\$debug);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --in not found.\n\n", -exitval => 2, -verbose => 1)  if (! $infile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --ref not found.\n\n", -exitval => 2, -verbose => 1)  if (! $ref);

## Check for MAFFT
my $MAFFT = `which mafft`;
unless ($MAFFT =~ m/mafft/) { die "\n Error! You need to make sure MAFFT is installed and included in your PATH before you can run this program...\n\n"; }

## Check that the input is peptide and NOT nucleotide
check_molecule($infile);

## Create the temporary working directory.
my @chars = ("A".."Z", "a".."z");
my $rand_string;
$rand_string .= $chars[rand @chars] for 1..8;
my $working_dir = "./762_caller_tmp_" . $rand_string;
print `mkdir -p $working_dir`;

###########
## Begin ##
###########
my $date = `date`;
print " 762-Caller Version: $version\n Beginning: $date";

## Read the references file into the scalar "References"
## I know, it seems crude, but it'll get the job done.
my $References;
my $nrefs=0;
open(IN,"<$ref") || die "\n Cannot open the file: $ref\n";
while(<IN>) {
    $References = $References . $_;
    if ($_ =~ m/^>/) { $nrefs++; }
}
close(IN);

## Read the FASTA into memory...
my %Fasta = ReadFasta($infile);
my @Order = FastaOrder($infile);



open(OUT,">$outfile") || die "\n Cannot write to the output file: $outfile\n";
print OUT "sequence\t762_residue\tpos-762\tpos-547\tpos-850\tcomplete?\n";
foreach my $header (@Order) {
    unless (exists $Fasta{$header}) { die "\n Missing sequence: $header\n\n"; }
    my $sequence = $Fasta{$header};
    open(TMP,">$working_dir/query_and_ref_tmp.fa") || die "Can't open the temporary FSA: $working_dir/query_and_ref_tmp.fa\n\n";
    print TMP $References . ">" . $header . "\n" . $sequence . "\n";
    close(TMP);
    # print `mafft --retree 2 --inputorder $working_dir/query_and_ref_tmp.fa > $working_dir/aln.fa 2> $working_dir/std.err`; ## Default MAFFT run, not local or global
    print `mafft --localpair  --maxiterate 16 --inputorder $working_dir/query_and_ref_tmp.fa > $working_dir/aln.fa 2> $working_dir/std.err`;
    flaten_fasta("$working_dir/aln.fa", "$working_dir/aln.flt.fa");
    my @results = parse_762("$working_dir/aln.flt.fa", $nrefs);
    print OUT join ("\t", @results) . "\n";
}
close(OUT);

## Clean up time!
unless ($debug) {
    if (-d $working_dir) {
	print `rm -rf $working_dir`;
    }
}
$date = `date`;
print " Done: $date";

## Subroutines
sub ReadFasta
{
    my $file = $_[0];
    $/=">";
    my %NTs;
    my @FILE;
    open(FASTA,"<$file") or die "\n\n\n Nada $file\n\n\n";
    @FILE=<FASTA>;
    close(FASTA);
    shift(@FILE);
    foreach my $orf (@FILE)
    {       my @Lines = split(/\n/,$orf);
	    my $name = $Lines[0];
	    my $seq = "";
	    foreach my $i (1..$#Lines)
	    {       $seq .= $Lines[$i]; }
	    $seq =~ s/>//;
	    $NTs{$name} = $seq;
    }
    $/="\n"; # reset input break character
    return %NTs;
}
sub FastaOrder
{
    my $file = $_[0];
    my @a;
    open(IN,"<$file") || die "\n Cannot open the file: $file\n";
    while(<IN>) {
	chomp;
	if ($_ =~ m/^>/) {
	    my $h = $_;
	    $h =~ s/^>//;
	    push(@a, $h);
	}
    }
    close(IN);
    return @a;
}
sub parse_762
{
    my $in = $_[0];
    my $query = $_[1] * 2;
    my @results;
    my $l=0;
    my ($header, $residue, $pos762, $pos547, $pos850, $complete) = ("","none",0,0,0,"no");
    my ($r547, $r762, $r850);
    my $pos=0;
    my $qlen;
    open(IN,"<$in") || die "\n Cannot open the file: $in\n";
    while(<IN>) {
	chomp;
	if ($l == 1) { # if on the E. coli sequence line...
	    my @A = split(//, $_);
	    for(my $i=0;$i<scalar(@A); $i++) {
		unless ($A[$i] eq "-") {
		    $pos++;
		    if ($pos == 547) { $r547 = $i; }
		    if ($pos == 762) { $r762 = $i; }
		    if ($pos == 850) { $r850 = $i; }
		}
	    }
	}
	elsif ($l == $query) {
	    $header = $_;
	    $header =~ s/^>//;
	}
	elsif ($l == $query+1) { # if on the query sequence line...
	    my @A = split(//, $_);
	    my $local_pos = 0;
	    for(my $i=0;$i<=$r850; $i++) {
		unless ($A[$i] eq "-") {
		    if ($i<=$r547) { $pos547++; }
		    if ($i<=$r762) { $pos762++; }
		    $pos850++;
		}
	    }
	    my $seq = $_;
	    $seq =~ s/-//g;
	    $qlen = length($seq);
	    unless ($pos762 == 0) { $residue = $A[$r762]; }
	    if ($pos547 != 0 && $pos850 != $qlen) { $complete = "yes"; }
	}
	$l++;
    }
    close(IN);
    @results = ($header, $residue, $pos762, $pos547, $pos850, $complete);
    return (@results);
}
sub flaten_fasta
{
    my $i = $_[0];
    my $o = $_[1];
    my $l=0;
    open(IN,"<$i")  || die "\n Cannot open the file: $i\n";
    open(OUTPUT,">$o") || die "\n Cannot write to: $o\n";
    while(<IN>) {
	chomp;
	if ($l==0) {
	    print OUTPUT $_ . "\n";
	}
	elsif ($_ =~ m/^>/) {
	    print OUTPUT "\n" . $_ . "\n";
	}
	else {
	    print OUTPUT $_;
	}
	$l++;
    }
    close(IN);
    print OUTPUT "\n";
    close(OUTPUT);
}
sub check_molecule
{
    my $i = $_[0];
    my $bases=0;
    my $ATGC=0;
    open(CHECK, "<$i") || die "\n Error: Cannot open the infile: $i\n";
    while(<CHECK>) {
	chomp;
	unless ($_ =~ m/^>/) {
	    $bases += length($_);
	    my $seq = $_;
	    $ATGC += $seq =~ tr/ATGCatgc/ATGCATGC/;
	}
    }
    close(CHECK);
    my $fraction = $ATGC / $bases;
    if ($fraction > 0.5) {
	$fraction *= 100;
	die "\n\n ERROR: Looks like you have submitted a nucleotide FASTA file and not a peptide FASTA file. Remember that the 762_caller needs peptide sequences to make the 762 call. I'm assuming this is a nucleotide FASTA file becasue $fraction % of the bases in this file are either A, T, G, or C...\n\n In the event that this is not correct email Dan (dnasko 'at' udel 'dot' edu) and let him know\n\n";
    }
}

exit 0;
