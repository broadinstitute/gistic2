#!/usr/bin/perl
# extract NCBI to ENSEMBL ID table from UCSC Homo_sapiens.gene_info
use strict;
use warnings;
print STDOUT "locus_id\tensembl_id\n";
while(<>) {
    if (!m/^#.*/) {
	my @fields = split("\t");
	my $locus_id = $fields[1];
	my $dbXrefs = $fields[5];
	for ($dbXrefs) {
	    if (m/Ensembl:([A-Z0-9]+)/) {
		my $ensembl_id = $1;
		print STDOUT "$locus_id\t$ensembl_id\n";
#!	    } else {
#!		print STDERR "<< $dbXrefs >>\n";
	    }
	}
    }
}
