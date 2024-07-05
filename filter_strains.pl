#!/usr/bin/env perl 

## filtering script for all strains with at least 1 read mapped to them
## key features: (1) all genes with too many mismatches are discarded (see cutoffs below) 
## (2) all genes with too many human reads are discarded (see cutoffs below)
## (3) all blacklisted genes are discarded 
## after this, any strain that has at least 1 detected rRNA and 1 protein gene are retained 

use strict; 
use warnings; 

if (scalar @ARGV != 3) {
  print STDERR "Usage: ./annotate_mapped_reads2.pl <accession_summary_tsv> <annotated_fcounts_tsv> <gene_blacklist>\n";
  exit 1;
}

my $acc_summary = shift @ARGV; 
my $ann_fcounts = shift @ARGV; 
my $gene_blacklist = shift @ARGV; 

## change cutoffs here; for now, 3 mm per ~90 bp read for rRNA and 2 mm per same for protein is what we used 
my $mismatch_rrna_cutoff = 3.4; 
my $mismatch_protein_cutoff = 2.3; 
my $human_fraction = 0.1; 

open ACCSUM,"<",$acc_summary or die "$!"; 
open FCOUNTS,"<",$ann_fcounts or die "$!"; 
open BLIST,"<",$gene_blacklist or die "$!"; 
open FSUM,">","filtered.summary.tsv" or die "$!"; 
open FGENE,">","filtered.gene.list" or die "$!"; 

my $BLG = {}; ## blacklisted genes
my $A2G = {}; ## gene properties, linked to assembly the gene belongs to

while (<BLIST>) {
  chomp; 
  $BLG->{$_} = 1; 
} 

while (<FCOUNTS>) {
  my @t = split /\t/;
  my $gene = $t[0]; 
  my $acc = $t[11]; 
  if (defined $A2G->{$acc}->{$gene}) { 
    print STDERR "WARNING: multiple lines per gene in annotated.fcounts.tsv! This should not happen; please investigate."; 
  } else {  
    $A2G->{$acc}->{$gene}->{type} = $t[8]; 
    $A2G->{$acc}->{$gene}->{mismatch} = $t[9]; 
    $A2G->{$acc}->{$gene}->{human} = $t[10];
	$A2G->{$acc}->{$gene}->{blacklisted} = 1 if (defined $BLG->{$gene}); 
  }
} 

while (<ACCSUM>) {
  my $acc = (split /\t/)[0];
  my $name = (split /\t/)[1];
  my ($prot,$rrna,$mis,$hum,$blg) = ('0') x 5; 
  foreach my $gene (keys %{$A2G->{$acc}}) { 
    if (defined $A2G->{$acc}->{$gene}->{blacklisted}) { 
	  $blg++; 
	  next; 
	} elsif ($A2G->{$acc}->{$gene}->{human} > $human_fraction) { 
	  $hum++; 
	  next; 
	} elsif ($A2G->{$acc}->{$gene}->{type} eq "rRNA" && $A2G->{$acc}->{$gene}->{mismatch} > $mismatch_rrna_cutoff) { 
	  $mis++;
	  next; 
	} elsif ($A2G->{$acc}->{$gene}->{type} ne "rRNA" && $A2G->{$acc}->{$gene}->{mismatch} > $mismatch_protein_cutoff) {
	  $mis++;
	  next; 
	} elsif ($A2G->{$acc}->{$gene}->{type} eq "rRNA") { 
	  $rrna++;
	  $A2G->{$acc}->{$gene}->{filtered} = 1;  
	  next; 
	} else { 
	  $prot++;
	  $A2G->{$acc}->{$gene}->{filtered} = 1;  
	}
  }
  ## finally, filter the accession table!
  if ($prot + $rrna >= 3) { 
    print STDOUT "Strain $acc ($name) is RETAINED: $rrna rRNA and $prot non-rRNA genes detected; not considered: $blg blacklisted, $hum potentially human-derived, and $mis with too many mismatches.\n"; 
	print FSUM;
	foreach my $gene (keys %{$A2G->{$acc}}) { 
	  print FGENE "$gene\n" if (defined $A2G->{$acc}->{$gene}->{filtered});
	} 
  } else {  
    print STDOUT "Strain $acc ($name) is REMOVED: $rrna rRNA and $prot non-rRNA genes detected; not considered: $blg blacklisted, $hum potentially human-derived, and $mis with too many mismatches.\n"; 
  } 
} 

close ACCSUM; 
close FCOUNTS; 
close BLIST; 
close FSUM;
close FGENE;
