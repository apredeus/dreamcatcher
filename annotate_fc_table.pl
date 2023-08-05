#!/usr/bin/env perl

## V2: using seq2taxid.map file from KrakenUniq
## filter non-zero genes, annotate w type and species, and make a summarised table on species level 

use strict; 
use warnings; 

if (scalar @ARGV != 3) { 
  print STDERR "Usage: ./annotate_fc_table.pl <fcounts_tsv> <gene_type_tsv> <seq2taxid.map>\n"; 
  exit 1; 
} 

my $fcounts_tsv = shift @ARGV; 
my $gene_type = shift @ARGV; 
my $seq2taxid = shift @ARGV; 

open TSV,"<",$fcounts_tsv or die "$!"; 
open TYPES,"<",$gene_type or die "$!"; 
open SEQ2TAXID,"<",$seq2taxid or die "$!"; 

open ANNTSV,">","annotated.fcounts.tsv"; 
open SUMMARY,">","accession.summary.tsv"; 

my $G2T = {}; 
my $C2A = {};
my $SUM = {}; 

print STDERR "Reading gene-to-gene type relationship table..\n"; 
while (<TYPES>) {
  chomp; 
  my @t = split /\t/; 
  $G2T->{$t[0]} = $t[1]; 
}

print STDERR "Reading seq2taxid table..\n";
while (<SEQ2TAXID>) {
  chomp;
  my @t = split /\t/;
	$t[2] =~ m/(GCF_\d+\.\d+) (.*?)$/;
  $C2A->{$t[0]}->{taxid} = $t[1]; 
  $C2A->{$t[0]}->{acc} = $1; 
  $C2A->{$t[0]}->{name} = $2; 
} 

print STDERR "Annotating the featureCounts output..\n"; 
while (<TSV>) { 
  chomp;
  next if (m/^#/ || m/^Geneid/); 
  my @t = split /\t/; 
  my $gene = $t[0]; 
  my $chr = $t[1];
  my $count = $t[6]; 
  if ($chr =~ m/(.*?);/) {
    print STDERR "WARNING: multiple chromosomes ($chr) for gene $gene! Will split them and use the first ID.\n"; 
    $chr = $1;
  }

  my ($acc,$name,$taxid,$type) = ("NONE") x 4; 

  if (defined $C2A->{$chr}) { 
    $acc = $C2A->{$chr}->{acc}; 
    $taxid = $C2A->{$chr}->{taxid}; 
    $name = $C2A->{$chr}->{name}; 
  } else { 
    print STDERR "WARNING: no accession/taxid/species defined for chromosome $chr! Reporting NONE.\n"; 
  } 
  if (defined $G2T->{$gene}) { 
    $type = $G2T->{$gene};
  } else { 
    print STDERR "WARNING: no type defined for gene $gene! Reporting NONE.\n"; 
  } 
  ## only print if read count > 0 
  if ($count > 0) {
    print ANNTSV "$_\t$type\t$acc\t$name\t$taxid\n";
    if (defined $SUM->{$acc}->{taxid}) {
      if ($type eq "rRNA") {
        $SUM->{$acc}->{rRNA_genes} = $SUM->{$acc}->{rRNA_genes} + 1;
        $SUM->{$acc}->{rRNA_count} = $SUM->{$acc}->{rRNA_count} + $count;
      } else { 
        $SUM->{$acc}->{protein_genes} = $SUM->{$acc}->{protein_genes} + 1;
        $SUM->{$acc}->{protein_count} = $SUM->{$acc}->{protein_count} + $count;
      } 
    } else { 
      $SUM->{$acc}->{name} = $name; 
      $SUM->{$acc}->{taxid} = $taxid; 
      if ($type eq "rRNA") { 
        $SUM->{$acc}->{rRNA_genes} = 1;
        $SUM->{$acc}->{protein_genes} = 0;
        $SUM->{$acc}->{rRNA_count} = $count;
        $SUM->{$acc}->{protein_count} = 0;
      } else { 
        $SUM->{$acc}->{protein_genes} = 1;
        $SUM->{$acc}->{rRNA_genes} = 0;
        $SUM->{$acc}->{protein_count} = $count;
        $SUM->{$acc}->{rRNA_count} = 0;
      }
    } 
  } 
} 

foreach my $acc (keys %{$SUM}) {
  my $taxid = $SUM->{$acc}->{taxid}; 
  my $name = $SUM->{$acc}->{name}; 
  my $rg = $SUM->{$acc}->{rRNA_genes}; 
  my $rc = $SUM->{$acc}->{rRNA_count}; 
  my $pg = $SUM->{$acc}->{protein_genes}; 
  my $pc = $SUM->{$acc}->{protein_count}; 

  print SUMMARY "$acc\t$taxid\t$name\t$rg\t$rc\t$pg\t$pc\n"; 
} 

close TSV; 
close TYPES; 
close SEQ2TAXID; 
close ANNTSV; 
close SUMMARY;
