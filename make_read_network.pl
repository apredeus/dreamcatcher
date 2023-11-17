#!/usr/bin/env perl 

## v2 - use GCF IDS for col1 & 2

use strict; 
use warnings;

if (scalar @ARGV != 2) { 
  print STDERR "Usage: ./make_read_network.pl <read_to_chr_filtered_uniq> <seq2taxid.map>\n"; 
  exit 1; 
}

my $read_to_chr = shift @ARGV;
my $seq2taxid = shift @ARGV; 

open READS,"<",$read_to_chr or die "$!"; 
open SEQ2TAXID,"<",$seq2taxid or die "$!"; 

my $C2A = {};    ## chrom-to-accession 
my $N = {};      ## nodes
my $E = {};      ## edges 

print STDERR "Reading seq2taxid table..\n";
while (<SEQ2TAXID>) {
	m/(GCF_\d+\.\d+)/;
	my $gcf = $1; 
	my $chr = (split /\t/)[0]; 
  $C2A->{$chr} = $gcf; 
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

  ## this is the simple trick to avoid putting everything into memory.
  if ($curr_read eq $prev_read) {
    push @read_chrs,$curr_chr;
  } elsif ($prev_read eq "") {
    $prev_read = $curr_read; 
    push @read_chrs,$curr_chr; 
  } else {
    my @read_gcfs;
    foreach my $chr (@read_chrs) {
      push @read_gcfs,$C2A->{$chr};
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
    @read_chrs = ();
    push @read_chrs,$curr_chr;
  } 

  $prev_read = $curr_read; 
}

## need extra thing for the last line of the file
my @read_gcfs;
foreach my $chr (@read_chrs) {
  push @read_gcfs,$C2A->{$chr};
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

close SEQ2TAXID;
close READS;

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
