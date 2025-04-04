#!/usr/bin/env perl

## V3: using seq2taxid.map file from KrakenUniq
## filter non-zero genes, annotate w type and species, and make a summarised table on species level
## v3 also uses per-read mismatch table to calculate mean mismatch per rRNA and protein coding genes

use strict;
use warnings;

if (scalar @ARGV != 5) { 
  print STDERR "Usage: ./make_summary_tables.pl <fcounts_tsv> <gene_type_tsv> <annotated_reads_tsv> <seq2taxid.map> <name_tag>\n"; 
  exit 1; 
} 

my $fcounts_tsv = shift @ARGV; 
my $gene_type = shift @ARGV; 
my $ann_reads = shift @ARGV; 
my $seq2taxid = shift @ARGV; 
my $name_tag = shift @ARGV; 

open FCOUNTS,"<",$fcounts_tsv or die "$!"; 
open TYPES,"<",$gene_type or die "$!"; 
open READS,"<",$ann_reads or die "$!"; 
open SEQ2TAXID,"<",$seq2taxid or die "$!"; 

open ANNFC,">",$name_tag.".annotated_fcounts.tsv"; 
open SUMMARY,">",$name_tag.".summary.tsv"; 

my $GENE = {}; 
my $CHR2ACC = {};
my $SUM = {}; 

print STDOUT "\tmake_summary_tables.pl: reading gene/gene_type relationship table..\n"; 
while (<TYPES>) {
  chomp; 
  my @t = split /\t/; 
  $GENE->{$t[0]}->{type} = $t[1];
}

print STDOUT "\tmake_summary_tables.pl: reading seq2taxid table..\n";
while (<SEQ2TAXID>) {
  chomp;
  my @t = split /\t/;
  $t[2] =~ m/(GCF_\d+\.\d+) (.*?)$/;
  $CHR2ACC->{$t[0]}->{taxid} = $t[1]; 
  $CHR2ACC->{$t[0]}->{acc} = $1; 
  $CHR2ACC->{$t[0]}->{name} = $2; 
} 

print STDOUT "\tmake_summary_tables.pl: reading read table, counting mismatches per gene..\n"; 
while (<READS>) { 
  chomp; 
  my @t = split /\t/;
  if ($t[1] eq "-") { 
    next; 
  } else {
    my @q = split /,/,$t[1]; ## sometimes one read is assigned to several genes
	## this part summarizes counts, human counts, and all mismatches
	foreach my $gene (@q) { 
      if (defined $GENE->{$gene}->{raw}) { 
	    $GENE->{$gene}->{raw} += 1; 
		$GENE->{$gene}->{mm} += $t[2]; 
		if ($t[3] ne "-" && $t[3] <= $t[2]) {
		  $GENE->{$gene}->{human} += 1;
		} 
	  } else { 
	    $GENE->{$gene}->{raw} = 1; 
		$GENE->{$gene}->{mm} = $t[2]; 
		if ($t[3] ne "-" && $t[3] <= $t[2]) {
		  $GENE->{$gene}->{human} = 1;
		} else { 
		  $GENE->{$gene}->{human} = 0;
		} 
	  }
    }
  } 
} 
	  
print STDOUT "\tmake_summary_tables.pl: annotating the featureCounts output..\n"; 
while (<FCOUNTS>) { 
  chomp;
  next if (m/^#/ || m/^Geneid/); 
  my @t = split /\t/; 
  my $gene = $t[0]; 
  my $chr = $t[1];
  my $fc_count = $t[6]; 
  if ($chr =~ m/(.*?);/) {
    print STDERR "\tmake_summary_tables.pl: WARNING: multiple chromosomes ($chr) for gene $gene! Will split them and use the first ID.\n"; 
    $chr = $1;
  }

  my ($acc,$name,$taxid,$type,$raw_count,$mean_mm,$human_frac) = ("NONE") x 7; 

  if (defined $CHR2ACC->{$chr}) { 
    $acc = $CHR2ACC->{$chr}->{acc}; 
    $taxid = $CHR2ACC->{$chr}->{taxid}; 
    $name = $CHR2ACC->{$chr}->{name}; 
  } else { 
    print STDERR "\tmake_summary_tables.pl: WARNING: no accession/taxid/species defined for chromosome $chr! Reporting NONE.\n"; 
  } 
  if (defined $GENE->{$gene}) { 
    $type = $GENE->{$gene}->{type};
	$raw_count = (defined $GENE->{$gene}->{raw}) ? $GENE->{$gene}->{raw} : 0;
  } else { 
    print STDERR "\tmake_summary_tables.pl: WARNING: no type defined for gene $gene! Reporting NONE.\n"; 
  } 
  ## only print if read count > 0 
  if ($raw_count > 0) {
    ## only calculate further stats for genes that are detected
    $mean_mm = $GENE->{$gene}->{mm}/$GENE->{$gene}->{raw};
	$human_frac = $GENE->{$gene}->{human}/$GENE->{$gene}->{raw};
    printf ANNFC "%s\t%d\t%s\t%.3f\t%.3f\t%s\t%s\t%d\n",$_,$raw_count,$type,$mean_mm,$human_frac,$acc,$name,$taxid;
    
	if (defined $SUM->{$acc}->{taxid}) {
      if ($type eq "rRNA") {
        $SUM->{$acc}->{rRNA_genes} += 1;
        $SUM->{$acc}->{rRNA_count} += $fc_count;
        $SUM->{$acc}->{rRNA_raw} += $raw_count;
        $SUM->{$acc}->{rRNA_mm} += $GENE->{$gene}->{mm};
        $SUM->{$acc}->{rRNA_human} += $GENE->{$gene}->{human};
      } else { 
        $SUM->{$acc}->{protein_genes} += 1;
        $SUM->{$acc}->{protein_count} += $fc_count;
        $SUM->{$acc}->{protein_raw} += $raw_count;
        $SUM->{$acc}->{protein_mm} += $GENE->{$gene}->{mm};
        $SUM->{$acc}->{protein_human} += $GENE->{$gene}->{human};
      } 
    } else { 
      $SUM->{$acc}->{name} = $name; 
      $SUM->{$acc}->{taxid} = $taxid;
      if ($type eq "rRNA") {
        $SUM->{$acc}->{rRNA_genes} = 1;
        $SUM->{$acc}->{rRNA_count} = $fc_count;
        $SUM->{$acc}->{rRNA_raw} = $raw_count;
        $SUM->{$acc}->{rRNA_mm} = $GENE->{$gene}->{mm};
        $SUM->{$acc}->{rRNA_human} = $GENE->{$gene}->{human};
      } else { 
        $SUM->{$acc}->{protein_genes} = 1;
        $SUM->{$acc}->{protein_count} = $fc_count;
        $SUM->{$acc}->{protein_raw} = $raw_count;
        $SUM->{$acc}->{protein_mm} = $GENE->{$gene}->{mm};
        $SUM->{$acc}->{protein_human} = $GENE->{$gene}->{human};
      } 
    } 
  } 
} 

foreach my $acc (keys %{$SUM}) {
  my $taxid = $SUM->{$acc}->{taxid}; 
  my $name = $SUM->{$acc}->{name};

  my ($rg,$rc,$rr,$rm,$rh,$pg,$pc,$pr,$pm,$ph) = ("0") x 10;
  if (defined $SUM->{$acc}->{rRNA_genes}) { 
    $rg = $SUM->{$acc}->{rRNA_genes}; 
    $rc = $SUM->{$acc}->{rRNA_count}; 
    $rr = $SUM->{$acc}->{rRNA_raw}; 
    $rm = $SUM->{$acc}->{rRNA_mm}/$SUM->{$acc}->{rRNA_raw}; 
    $rh = $SUM->{$acc}->{rRNA_human}/$SUM->{$acc}->{rRNA_raw}; 
  }
  if (defined $SUM->{$acc}->{protein_genes}) { 
    $pg = $SUM->{$acc}->{protein_genes}; 
    $pc = $SUM->{$acc}->{protein_count}; 
    $pr = $SUM->{$acc}->{protein_raw}; 
    $pm = $SUM->{$acc}->{protein_mm}/$SUM->{$acc}->{protein_raw}; 
    $ph = $SUM->{$acc}->{protein_human}/$SUM->{$acc}->{protein_raw}; 
  } 

  printf SUMMARY "%s\t%s\t%d\t%d\t%.6f\t%d\t%.3f\t%.3f\t%d\t%.6f\t%d\t%.3f\t%.3f\n",$acc,$name,$taxid,$rg,$rc,$rr,$rm,$rh,$pg,$pc,$pr,$pm,$ph; 
} 

close FCOUNTS; 
close TYPES; 
close READS;
close SEQ2TAXID;

close ANNFC; 
close SUMMARY;
