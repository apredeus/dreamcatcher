#!/usr/bin/env perl 

## we miss pseudogenes here. 99% sure it does not matter. 

use strict; 
use warnings; 

if (scalar @ARGV != 1) { 
    print STDERR "USAGE: ./extract_gene_biotypes.pl <detected_strains.gff>\n"; 
    exit 1; 
}

my $gff = shift @ARGV; 
open GFF,"<",$gff or die "$!"; 

while (<GFF>) {
    if (m/\tgene\t.*gene_biotype=(.*?);/) {
        my $type = $1;
        my $lt = ""; 
        if (m/;locus_tag=(.*?);/) {
            $lt = $1;
        } elsif (m/;locus_tag=(.*?)$/) {
            $lt = $1;
        }
        print "$lt\t$type\n";
    } 
} 

close GFF; 
