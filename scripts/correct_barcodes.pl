#!/usr/bin/env perl 

use strict;
use warnings; 

if (scalar @ARGV != 3 ) { 
	print STDERR "Usage: ./correct_barcodes.pl <se_umi_filt.fastq> <barcode_correction.tsv> <output fastq>\n"; 
	exit 1; 
} 

my $fastq = shift @ARGV; 
my $bc_tsv = shift @ARGV; 
my $out = shift @ARGV; 

open FASTQ,"<",$fastq or die "$!"; 
open TSV,"<",$bc_tsv or die "$!"; 
open OUT,">",$out or die "$!"; 

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
			print OUT; 
			$exact++; 
		} elsif (! defined $BC->{$bc}) { 
			print OUT; 
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
			print OUT; 
			$corrected++;
		}
	} else { 
		print OUT; 
	}
} 

printf STDOUT "\t\tcorrect_barcodes.pl: processed fastq file with %d total reads; %d matching a whitelist, %d corrected, and %d unable to correct!\n",$exact+$corrected+$unmatched,$exact,$corrected,$unmatched; 

close FASTQ; 
close TSV;
close OUT; 
