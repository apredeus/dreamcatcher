#!/usr/bin/env perl 

## remove genes that are to be removed, and add GCF strain ID to each gene 

use strict; 
use warnings; 

if (scalar @ARGV != 3) {
	print STDERR "Usage: ./filter_top_gff.pl <unfiltered_top_gff> <genes_to_remove> <seq2taxid.map>\n";
	exit 1;
}

my $top_gff = shift @ARGV; 
my $genes = shift @ARGV; 
my $seq2taxid = shift @ARGV; 

open GFF,"<",$top_gff or die "$!"; 
open GLIST,"<",$genes or die "$!"; 
open SEQ2TAXID,"<",$seq2taxid or die "$!"; 

my $BLG = {};
my $CHR2ACC = {};
my $removed = 0; 
my $kept = 0; 

while (<SEQ2TAXID>) {
    chomp;
    my @t = split /\t/;
    $t[2] =~ m/(GCF_\d+\.\d+) (.*?)$/;
    $CHR2ACC->{$t[0]}->{taxid} = $t[1];
    $CHR2ACC->{$t[0]}->{acc} = $1;
    $CHR2ACC->{$t[0]}->{name} = $2;
}

while (<GLIST>) { 
    chomp;
    $BLG->{$_} = 1; 
} 

while (<GFF>) { 
    my @t = split /\t/; 
    if (m/\tgene\t.*gene_biotype=(.*?);/) {
        my $type = $1;
        my $lt = "";
        my $chr = $t[0];
        my $begin = $t[3]; 
        my $end = $t[4]; 
        my $strand = $t[6]; 
        if ($t[8] =~ m/;locus_tag=(.*?);/) {
            $lt = $1;
        } elsif ($t[8] =~ m/;locus_tag=(.*?)$/) {
            $lt = $1;
        }
        my $gcf = $CHR2ACC->{$chr}->{acc};
        
        ## don't print if gene was removed 
        if (defined $BLG->{$lt}) { 
            $removed++; 
        } else { 
            $kept++; 
            print "$chr\tRefSeq\tgene\t$begin\t$end\t.\t$strand\t.\tlocus_tag=$lt;strain_id=$gcf;gene_biotype=$type;\n";
        }
    }
}

#print STDERR "\tfilter_top_gff.pl: kept $kept genes, removed $removed genes..\n"; 

close GLIST;
close SEQ2TAXID; 
close GFF; 
