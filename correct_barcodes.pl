#!/usr/bin/env perl 

use strict;
use warnings; 
use Data::Dumper; 

if (scalar @ARGV != 2 ) { 
	print STDERR "Usage: ./correct_barcodes.pl <se_umi_filt.fastq> <barcode_correction.tsv>\n"; 
	exit 1; 
} 

my $fastq = shift @ARGV; 
my $bc_tsv = shift @ARGV; 

open FASTQ,"<",$fastq or die "$!"; 
open TSV,"<",$bc_tsv or die "$!"; 

my $FQ = {}; 
my $BC = {}; 

while (<TSV>) {
	if (m/\t0$/) {
		my $bc = (split /\t/)[0]; 
		$BC->{$bc} = "exact"; 
	} elsif (m/\t1$/) { 
		my $bc = (split /\t/)[0];
		my $corr = (split /\t/)[1];
		if (defined $BC->{$bc} && $BC->{$bc} eq "exact") { 
			next; 
		} else { 
			push @{$BC->{$bc}}, $corr;
		} 
	} 
}

while (<FASTQ>) { 
	if (m/^@.*_([A-Z]+)_([A-Z]+)/) { 
		my $bc = $1;
		if (defined $BC->{$bc} && $BC->{$bc} eq "exact") { 
			$FQ->{$bc} += 1; 
		} 
	} 
} 
close FASTQ; 

open FASTQ,"<",$fastq or die "$!"; 
my ($exact,$corrected,$unmatched) = (0) x 3; 

while (<FASTQ>) {
	if (m/^@.*_([A-Z]+)_([A-Z]+)/) {
		my $bc = $1; 
		if (defined $BC->{$bc} && $BC->{$bc} eq "exact") { 
			print; 
			$exact++; 
		} elsif (! defined $BC->{$bc}) { 
			print; 
			$unmatched++; 
		} else { 
			my $corr = @{$BC->{$bc}}[0];

			## if there are multiple barcodes this one can be corrected to
			## find one that occurs in most reads, and use that one
			if (scalar @{$BC->{$bc}} > 1) { 
				my $hits = (defined $FQ->{$corr}) ? $FQ->{$corr} : 0; 
				foreach my $match (@{$BC->{$bc}}) { 
					if (defined $FQ->{$match} && $FQ->{$match} > $hits) { 
						$corr = $match; 
						$hits = $FQ->{$match}; 
					} 
				}
			} 
			s/_${bc}_/_${corr}_/;
			print; 
			$corrected++;
		}
	} else { 
		print; 
	}
} 

printf STDERR "Processed fastq file with %d total reads; %d matching a whitelist, %d corrected, and %d unable to correct!\n",$exact+$corrected+$unmatched,$exact,$corrected,$unmatched; 

close FASTQ; 
close TSV; 
