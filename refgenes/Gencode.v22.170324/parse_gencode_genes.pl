#!/usr/bin/perl
# filter gencode annotation.gtf to a standard TSV
use strict;
use warnings;
print STDOUT "chr\tsource\tstart\tend\tstrand\tgene_id\tgene_name\tgene_type\tgene_status\tlevel\n";
while(<>) {
    if (!m/^##.*/) {
	my @fields = split("\t");
	my $type = $fields[2];
	if ($type eq "gene") {
	    my $chr = $fields[0];
	    my $source = $fields[1];
	    my $start = $fields[3];
	    my $end = $fields[4];
	    my $strand = $fields[6];
	    my $attribute_str = $fields[8];
	    my @attribute_pairs = split(";",$attribute_str);
	    my %att = ();
	    foreach(@attribute_pairs) {
		if (m/(\S+)\s+"(\S+)"/) {
		    $att{$1} = $2;
		}
		else {
		    m/(\S+)\s+(\S+)/;
		    $att{$1} = $2;
		}
	    }
	    print STDOUT "$chr\t$source\t$start\t$end\t$strand\t$att{gene_id}\t$att{gene_name}\t$att{gene_type}\t$att{gene_status}\t$att{level}\n";
	}
    }
}
