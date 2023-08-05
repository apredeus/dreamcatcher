#!/usr/bin/env perl 

use strict; 
use warnings; 

if (scalar @ARGV != 3) { 
  print STDERR "Usage: ./make_read_network.pl <acc_to_species> <chr_to_acc> <read_to_chr>\n"; 
  exit 1; 
}

my $species_to_acc = shift @ARGV;
my $chr_to_acc = shift @ARGV; 
my $read_to_chr = shift @ARGV;

open SPECIES,"<",$species_to_acc or die "$!"; 
open CHR,"<",$chr_to_acc or die "$!"; 
open READS,"<",$read_to_chr or die "$!"; 

open NETWORK,">","nodes_and_edges.tsv"; 

my $C2A = {}; 
my $A2S = {};

my $N = {}; 
my $E = {}; 

print STDERR "Reading chromosome name to accession table..\n";
while (<CHR>) {
  chomp;
  my @t = split /\t/; 
  $C2A->{$t[0]} = $t[1]; 
} 
print STDERR "Reading accession to name/taxid table..\n"; 
while (<SPECIES>) {
  chomp; 
  my @t = split /\t/;
  ## there is an invisible empty column there
  $A2S->{$t[0]} = $t[3]; 
}

## take advantage of the read-sorted/line-uniq'd read_name-to-chr file 
my $prev_read = "";
my @read_chrs; 

print STDERR "Processing the reads..\n"; 
while (<READS>) { 
  chomp; 
  my @t = split /\t/;
  my $curr_read = $t[0]; 
  my $curr_chr = $t[1];

  ## this is the trick to avoid putting everything into memory.
  if ($curr_read eq $prev_read) {
    push @read_chrs,$curr_chr;
  } elsif ($prev_read eq "") {
    $prev_read = $curr_read; 
    push @read_chrs,$curr_chr; 
  } else {
    my @read_species;
    foreach my $chr (@read_chrs) {
      my $acc = $C2A->{$chr};
      my $species = $A2S->{$acc}; 
      push @read_species,$species;
    }

    my @uniq_species = uniq(@read_species);
    # printf STDERR "DEBUG: %s %d %s\n",$prev_read,scalar(@uniq_species),join(' ',@uniq_species); 
    foreach my $sp1 (@uniq_species) { 
      $N->{$sp1} += 1; 
      foreach my $sp2 (@uniq_species) { 
        $E->{$sp1}->{$sp2} += 1;
      } 
    }
    ## next time, flush (c) 
    @read_chrs = ();
    push @read_chrs,$curr_chr;
  } 

  $prev_read = $curr_read; 
}

## need extra thing for the last line of the file
my @read_species;
foreach my $chr (@read_chrs) {
  my $acc = $C2A->{$chr};
  my $species = $A2S->{$acc}; 
  push @read_species,$species;
}

## only uniq species (so sp1 will never be equal to sp2 below): 
my @uniq_species = uniq(@read_species); 
foreach my $sp1 (@uniq_species) { 
  $N->{$sp1} += 1; 
  foreach my $sp2 (@uniq_species) { 
    $E->{$sp1}->{$sp2} += 1;
  } 
}

printf STDERR "DEBUG: %s %d %s\n",$prev_read,scalar(@uniq_species),join(' ',@uniq_species); 

foreach my $sp1 (keys %{$E}) { 
  foreach my $sp2 (keys %{$E->{$sp1}}) {
    my $size1 = $N->{$sp1}; 
    my $size2 = $N->{$sp2};
    ## weighted overlap = overlap divided by smallest set size
    my $wt = ($size1 >= $size2) ? $E->{$sp1}->{$sp2}/$size2 : $E->{$sp1}->{$sp2}/$size1;
    ## 'gt' below only prints the values one side of the matrix diagonal 
    printf NETWORK "%s\t%s\t%d\t%d\t%d\t%.6f\n",$sp1,$sp2,$size1,$size2,$E->{$sp1}->{$sp2},$wt if ($wt >= 0.1 && $sp1 gt $sp2);
  } 
} 

close NETWORK; 
close SPECIES;
close CHR; 
close READS;

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
