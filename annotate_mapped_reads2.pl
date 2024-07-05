#!/usr/bin/env perl

## v2 of read annotation script calculates mismatch (NM) to mapped read length ratio
## this way it's more compareable across samples
## (this would be a bitch to implement for paired-end reads, wouldn't it)

use strict; 
use warnings; 

if (scalar @ARGV != 2) { 
  print STDERR "Usage: ./annotate_mapped_reads2.pl <fcounts_annotated_bam> <human_remap_bam>\n"; 
  exit 1; 
} 

my $fcounts_bam = shift @ARGV; 
my $human_bam = shift @ARGV; 

open FCBAM,"samtools view -\@4 $fcounts_bam |" or die "ERROR: failed to open combined bacterial bam $fcounts_bam using samtools!";
open HUMBAM,"samtools view -\@4 $human_bam |" or die "ERROR: failed to open human remap bam $human_bam file using samtools!";

my $H = {}; 

while (<HUMBAM>) {
  chomp; 
  my @t = split /\t/;
  m/NM:i:(\d+)/; 
  my $hum_mis = $1;
  my $mapped_length = 0;
  ## find the number of mapped bases from the CIGAR string
  my @maps = $t[5] =~ m/(\d+[A-Z])/g;
  foreach my $map (@maps) {
    $mapped_length += $1 if ($map =~ m/(\d+)M$/); 
  } 
  ## report mismatches per 100 mapped bp
  $H->{$t[0]} = sprintf '%.3f', 100*$hum_mis/$mapped_length; 
}

while (<FCBAM>) { 
  chomp;
  my @t = split /\t/; 
  m/NM:i:(\d+)/; 
  my $bac_mis = $1; 
  my $mapped_length = 0;
  ## find the number of mapped bases from the CIGAR string
  my @maps = $t[5] =~ m/(\d+[A-Z])/g;
  foreach my $map (@maps) {
    $mapped_length += $1 if ($map =~ m/(\d+)M$/); 
  } 
  my $bac_mmrate = sprintf '%.3f', 100*$bac_mis/$mapped_length; 
  my $rname = $t[0]; 
  my $gene = "-";
  my $hum_mmrate = (defined $H->{$rname}) ? $H->{$rname} : "-";
  $gene = $1 if (m/XS:Z:Assigned.*XT:Z:(.*?)$/); 

  print STDOUT "$rname\t$gene\t$bac_mmrate\t$hum_mmrate\n"; 
} 

close FCBAM; 
close HUMBAM;
