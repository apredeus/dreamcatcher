#!/usr/bin/env perl 

## takes extended taxDB from KrakenUniq 
## (with 'fake' taxids assigned to assemblies and sequences) 
## and selects only the relevant bacterial genera, species, and assemblies
## final result is a table of 5 columns: GCF, species name, species taxid, genus name, genus taxid

use strict; 
use warnings; 
use Data::Dumper; 

if (scalar @ARGV != 2) { 
	print "Usage: ./parse_taxdb.pl <extended_taxonomy_file> <Refseq GCF list>\n";
	exit 1;
}

my $taxdb = shift @ARGV; 
my $gcfs = shift @ARGV; 

open GCFS,"<",$gcfs or die "$!";
my $GCF = {};

print STDERR "Processing the list of allowed RefSeq IDs..\n";  
while (<GCFS>) { 
	chomp;
	## hash all relevant bacterial GCFs for a quick yes/no lookup
	$GCF->{$_} = 1; 
} 

## taxonomy is complicated, read this thing in child->parent pairs
## recent NCBI taxonomy fit into about 6 Gb
print STDERR "Processing the taxonomy file and making the pairwise hash..\n"; 
open TAXDB,"<",$taxdb or die "$!"; 
my $TAX = {};
my $lines = 0; 
while (<TAXDB>) {
	print STDERR "Processed $lines lines..\n" if ($lines % 200000 == 0);
	if (! m/\tassembly$/ && ! m/\tsequence$/) {
		chomp;
		my @t = split /\t/; 
		my $child_id = $t[0]; 
		my $parent_id = $t[1]; 
		my $child_name = $t[2]; 
		my $child_rank = $t[3]; 
		$TAX->{$child_id}->{rank} = $child_rank; 
		$TAX->{$child_id}->{parent} = $parent_id;
		$TAX->{$child_id}->{name} = $child_name;
	}
	$lines++; 
}
close TAXDB;
#print Dumper $TAX;

## now comes the tricky part! find an appropriate species and genus for all GCF
print STDERR "Finding the matching species/genus for each RefSeq ID provided..\n"; 
open TAXDB,"<",$taxdb or die "$!";
$lines = 0; 
while (<TAXDB>) { 
	if (m/\tassembly/ && m/\t(GCF_\d+\.\d+)/) {
	print STDERR "Processed $lines assembly IDs..\n" if ($lines % 2000 == 0);
		my $gcf = $1;
		if (defined $GCF->{$gcf}) {
			my $id = (split /\t/)[1];
			my ($species_id,$genus_id,$species_name,$genus_name) = (undef) x 4;
			my $cycles = 0; 
			while ((! defined $species_id || ! defined $genus_id ) && $id != '1' && $cycles < 100) { 
				$species_id = $id if ($TAX->{$id}->{rank} eq 'species'); 
				$genus_id = $id if ($TAX->{$id}->{rank} eq 'genus');
				$id = $TAX->{$id}->{parent};
				$cycles++; 
			}

  		if (defined $species_id) {
  			$species_name = $TAX->{$species_id}->{name};
  		} else { 
  			$species_id = "NONE"; 
  			$species_name = "NONE"; 
  		} 
  
  		if (defined $genus_id) { 
  			$genus_name = $TAX->{$genus_id}->{name}; 
  		} else { 
  			$genus_id = "NONE";
  			$genus_name = "NONE"; 
  		}

			print "$gcf\t$species_id\t$species_name\t$genus_id\t$genus_name\n"; 
		}
		$lines++;
	}
}

close TAXDB; 
