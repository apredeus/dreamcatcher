#!/usr/bin/env perl 

use strict; 
use warnings; 

if (scalar @ARGV != 2) { 
	print STDERR "USAGE: ./barcode_to_species.pl <gzipped_umitools_counts> <gcf_species_genus_tsv>\n"; 
	exit 1; 
} 

my $umi_out = shift @ARGV; 
my $gsa = shift @ARGV; 

open UMI,"zcat $umi_out|" or die "$!"; 
open GCF2TAXID,"<",$gsa or die "$!"; 

my $GSA = {}; 
my $BC = {}; 
my @sp = (); 

while (<GCF2TAXID>) { 
	my @t = split /\t/;
	$GSA->{$t[0]} = $t[2]; 
} 

while (<UMI>) { 
	next if (m/^gene\t/);
	chomp;
	my @t = split /\t/;
	my $gcf = $t[0]; 
	my $bc = $t[1]; 

	## if maps to > 1 strain, there will be a comma (there won't - since diff 
    ## strains have diff chromosomes, they will be present in BAM as separate lines) 
    $gcf =~ s/,.*//g;
	my $sp; 
	if (defined $GSA->{$gcf}) {
		$sp = $GSA->{$gcf};
        push @sp, $sp; 
	} else {
        print STDERR "barcode_to_species.pl: WARNING - no species found in the GSA file for $gcf! Setting to NONE.."; 
    } 

	if (! defined $BC->{$bc}->{$sp}) { 
		$BC->{$bc}->{$sp} = $t[2]; 
	}	else { 
		$BC->{$bc}->{$sp} += $t[2];
	} 
} 

## print header
my @uniq_sp = uniq(@sp);  
print "barcode"; 
foreach my $sp (@uniq_sp) { 
	print "\t$sp"; 
} 
print "\n"; 

## print per-barcode per-species UMI counts 
foreach my $bc (keys %{$BC}) {
	print $bc; 
	foreach my $sp (@uniq_sp) {
		my $count = (defined $BC->{$bc}->{$sp}) ? $BC->{$bc}->{$sp} : 0; 
		print "\t$count"; 
	} 
	print "\n"; 
} 

close UMI; 
close GCF2TAXID;

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
