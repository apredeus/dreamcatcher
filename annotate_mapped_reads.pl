#!/usr/bin/env perl

## V2: using seq2taxid.map file from KrakenUniq
## filter non-zero genes, annotate w type and species, and make a summarised table on species level 

use strict; 
use warnings; 

if (scalar @ARGV != 2) { 
  print STDERR "Usage: ./annotate_mapped_reads.pl <fcounts_annotated_bam> <human_read_mismatch_tsv>\n"; 
  exit 1; 
} 

my $fcounts_bam = shift @ARGV; 
my $human_mismatch = shift @ARGV; 

open BAM,"samtools view -\@4 $fcounts_bam |" or die "ERROR: failed to open the bam file using samtools!";
open HUMAN,"<",$human_mismatch or die "$!"; 

my $H = {}; 

while (<HUMAN>) {
  chomp; 
  my @t = split /\t/; 
  $H->{$t[0]} = $t[1]; 
}

while (<BAM>) { 
  chomp;
  my @t = split /\t/; 
  m/NM:i:(\d+)/; 
  my $bac_mis = $1; 
  my $rname = $t[0]; 
  my $gene = "-";
  my $human_mis = (defined $H->{$rname}) ? $H->{$rname} : "-";
  $gene = $1 if (m/XS:Z:Assigned.*XT:Z:(.*?)$/); 

  print STDOUT "$rname\t$gene\t$bac_mis\t$human_mis\n"; 
} 

close BAM; 
close HUMAN;
