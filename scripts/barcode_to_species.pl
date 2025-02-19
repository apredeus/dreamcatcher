#!/usr/bin/env perl 

use strict; 
use warnings; 

if (scalar @ARGV != 2) { 
	print STDERR "USAGE: ./barcode_to_species.pl <gzipped_umitools_counts> <prefix_to_species_tsv>\n"; 
	exit 1; 
} 

my $umi_out = shift @ARGV; 
my $prefix = shift @ARGV; 

open UMI,"zcat $umi_out|" or die "$!"; 
open PREFIX,"<",$prefix; 

my $PF = {}; 
my $BC = {}; 

while (<PREFIX>) { 
	chomp; 
	my @t = split /\t/;
	$PF->{$t[0]} = $t[1]; 
} 

while (<UMI>) { 
	next if (m/^gene\t/);
	chomp;
	my @t = split /\t/;
	my $pf = $t[0]; 
	my $bc = $t[1]; 

	$pf =~ s/,.*//g; ## reads split between 2 reads are reported by featureCounts like this: BI380_RS19515,BI380_RS19520
	$pf =~ s/\d+$//g;

	## the following is just a stupid hack to a stupid problem (multiple prefixes per 1 annotation, eg ECs_/ECs_R)
	my $sp; 
	if (defined $PF->{$pf}) {
		$sp = $PF->{$pf};
	} else { 
		foreach my $kpf (keys %{$PF}) { 
			$pf = $kpf if ($pf =~ m/^$kpf/); 
		}
		$sp = $PF->{$pf};
	} 
	if (! defined $BC->{$bc}->{$sp}) { 
		$BC->{$bc}->{$sp} = $t[2]; 
	}	else { 
		$BC->{$bc}->{$sp} += $t[2];
	} 
} 

## print header
my @sp = values %{$PF}; 
print "barcode"; 
foreach my $sp (@sp) { 
	print "\t$sp"; 
} 
print "\n"; 

## print per-barcode per-species UMI counts 
foreach my $bc (keys %{$BC}) {
	print $bc; 
	foreach my $sp (@sp) {
		my $count = (defined $BC->{$bc}->{$sp}) ? $BC->{$bc}->{$sp} : 0; 
		print "\t$count"; 
	} 
	print "\n"; 
} 

close UMI; 
close PREFIX; 
