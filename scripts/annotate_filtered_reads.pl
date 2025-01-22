#!/usr/bin/env perl

## take the filtered tsv file, add gene type & 
## also split the commas (multigene spanning reads) 

use strict; 
use warnings; 

if (scalar @ARGV != 2) { 
  print STDERR "Usage: ./annotate_filtered_reads.pl <annotated.fcounts.tsv> <read_to_gene_filtered_uniq.tsv>\n"; 
  exit 1; 
} 

my $ann_fcounts = shift @ARGV; 
my $filt_reads = shift @ARGV; 

open ANNFC,"<",$ann_fcounts or die "$!"; 
open READS,"<",$filt_reads or die "$!"; 

my $H = {}; 

while (<ANNFC>) {
  my @t = split /\t/;
  $H->{$t[0]}->{type} = $t[8]; 
  $H->{$t[0]}->{refseq} = $t[11]; 
}

while (<READS>) { 
  chomp;
  my $read = (split /\t/)[0];
  my $gene = (split /\t/)[1];
  if ($gene !~ m/,/) { 
    printf STDOUT "%s\t%s\t%s\t%s\n",$read,$gene,$H->{$gene}->{type},$H->{$gene}->{refseq};
  } else { 
    my @genes = split /,/, $gene; 
    foreach my $spg (@genes) { 
      printf STDOUT "%s\t%s\t%s\t%s\n",$read,$spg,$H->{$spg}->{type},$H->{$spg}->{refseq};
    }
  } 
} 

close ANNFC; 
close READS;
