#!/usr/bin/env perl 

## v3 - use GCF IDS, and only reads assigned to a gene

use strict; 
use warnings;

if (scalar @ARGV != 2) { 
  print STDERR "Usage: ./make_read_network2.pl <read_to_gene_filtered_uniq> <filtered_fcounts_tsv>\n"; 
  exit 1; 
}

my $read_to_gene = shift @ARGV;
my $filt_fc = shift @ARGV; 

open READS,"<",$read_to_gene or die "$!"; 
open FILTFC,"<",$filt_fc or die "$!"; 

my $G2A = {};    ## gene-to-accession 
my $N = {};      ## nodes
my $E = {};      ## edges 

#print STDOUT "make_read_network.pl: reading seq2taxid table..\n";
while (<FILTFC>) {
  my $gene = (split /\t/)[0]; 
  my $gcf = (split /\t/)[11]; 
  $G2A->{$gene} = $gcf; 
} 

## take advantage of the read-sorted/line-uniq'd read_name-to-gene file 
my $prev_read = "";
my @read_genes; 

#print STDOUT "make_read_network.pl: processing the reads..\n"; 
while (<READS>) { 
  chomp; 
  my @t = split /\t/;
  my $curr_read = $t[0]; 
  my $curr_gene;
  ## in rare cases, there would be multiple comma-separated genes the read mapped to, but only 1 would be in a filtered list
  if ($t[1] =~ m/,/) { 
    my @q = split /,/,$t[1];
    foreach my $gene (@q) {
      if (defined $G2A->{$gene}) {
	    $curr_gene = $gene;
	    last;
	  }
    } 
  } else { 
    $curr_gene = $t[1]; 
  } 

  ## this is the simple trick to avoid putting everything into memory.
  if ($curr_read eq $prev_read) {
    push @read_genes,$curr_gene;
  } elsif ($prev_read eq "") {
    $prev_read = $curr_read; 
    push @read_genes,$curr_gene;
  } else {
    my @read_gcfs;
    foreach my $gene (@read_genes) {
      push @read_gcfs,$G2A->{$gene} if (defined $G2A->{$gene});
    }

    my @uniq_gcfs = uniq(@read_gcfs);
    # printf STDERR "DEBUG: %s %d %s\n",$prev_read,scalar(@uniq_species),join(' ',@uniq_species); 
    foreach my $gcf1 (@uniq_gcfs) { 
      $N->{$gcf1} += 1; 
      foreach my $gcf2 (@uniq_gcfs) { 
        $E->{$gcf1}->{$gcf2} += 1;
      } 
    }
    ## next time, flush (c) 
    @read_genes = ();
    push @read_genes,$curr_gene;
  } 

  $prev_read = $curr_read; 
}

## need extra thing for the last line of the file
my @read_gcfs;
foreach my $gene (@read_genes) {
  push @read_gcfs,$G2A->{$gene} if (defined $G2A->{$gene});
}

## only uniq strains (so strain1 will never be equal to strain2 below): 
my @uniq_gcfs = uniq(@read_gcfs); 
foreach my $gcf1 (@uniq_gcfs) {
  $N->{$gcf1} += 1; 
  foreach my $gcf2 (@uniq_gcfs) { 
    $E->{$gcf1}->{$gcf2} += 1;
  } 
}

foreach my $gcf1 (keys %{$E}) { 
  foreach my $gcf2 (keys %{$E->{$gcf1}}) {
    my $size1 = $N->{$gcf1}; 
    my $size2 = $N->{$gcf2};
    ## weighted overlap = overlap divided by smallest set size
    my $wt = ($size1 >= $size2) ? $E->{$gcf1}->{$gcf2}/$size2 : $E->{$gcf1}->{$gcf2}/$size1;
    ## 'gt' below only prints the values one side of the matrix diagonal 
		## in this version of the script, we don't filter by $wt to see the distribution
    printf "%s\t%s\t%d\t%d\t%d\t%.6f\n",$gcf1,$gcf2,$size1,$size2,$E->{$gcf1}->{$gcf2},$wt if ($gcf1 ge $gcf2);
  } 
} 

close FILTFC;
close READS;

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
