![Alt text](https://github.com/dnasko/dna_pola_762_caller/blob/master/images/logo.png "762!")

Welcome to the viral DNA Polymerase A *E. coli* 762 caller.

This program was inspired by Schmidt et al. [Shotgun metagenomics indicates novel family A DNA polymerases predominate within marine virioplankton](http://www.nature.com/ismej/journal/v8/n1/full/ismej2013124a.html).

In short, the study concluded that when aligning viral DNA polymerase A proteins with *E. coli* DNA polymerase A, the amino acid present in a virus's DNA polymerase relative to *E. coli*'s 762nd amino acid can be indicative of viral lifestyle (i.e. lytic, lysogenic). The wildtype is for a phenylalanine (Phe) to be present at the 762 position, however Schmidt et al observed two distinct mutations that were also present in dsDNA viruses from the environment: leucine (Leu) and tyrosine (Tyr). Below is a table briefly desrcibing the implications of each mutation:

| Mutation      | Predicted Viral Lifestyle | 
|---------------|---------------------------|
| Phenylalanine | Lytic                     | 
| Leucine       | Potentially lysogenic     | 
| Tyrosine      | Lytic                     | 

It has been demonstrated that a leucine substitution at position 762 leads to a slower, but higher fidelity, DNA polymerase.

Usage
-----

> ./762_caller.pl --help

Given an input FASTA of viral DNA polymerase A peptides the 762_caller will perform iterative multiple sequence alignments (MSA) using [MAFFT](http://mafft.cbrc.jp/alignment/software/) against a set of reference sequences; cheif amung them *E. coli*.

The 762_caller will output a 6 column tab-delimmited text file with one row of results for each input sequence, each of the 6 fields are detailed below:

1. Sequence ID
2. 762 Position relative to this sequence
3. Residue at this sequence's 762
4. 547 Position relative to this sequence
5. 926 Position relative to this sequence
6. Does the sequence span the 547 - 926 trimming regions

Post-processing
---------------

Once you have your output from the 762_caller you can use trim_sequences.pl in the [scripts directory](https://github.com/dnasko/dna_pola_762_caller/tree/master/scripts) to trim your sequences to the 547 -> 926 region, as well as filter out any sequences that don't span this region end-to-end.

Acknowledgements
----------------

Support from the University of Delaware Center for Bioinformatics and Computational Biology Core Facility and use of the BIOMIX compute cluster was made possible through funding from Delaware INBRE (NIGMS GM103446) and the Delaware Biotechnology Institute.

DNA polymerase image credit niehs.nih.gov
