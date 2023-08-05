#!/usr/bin/env perl 

## this script will take Kraken2 and/or KrakenUniq report, extract the taxonomy IDs *for bacteria*, 
## and match the ID to a RefSeq genome ID from seq2taxid map file 
## in cases where one taxid matches many genomes, representative genome will be used (this requires a separate file)

use strict; 
use warnings; 
use Data::Dumper; 

if (scalar @ARGV != 3) { 
	print "USAGE: ./kraken_to_GCF.pl <KrakenUniq_report> <seq2taxid> <genbank_acc_file>\n"; 
	exit 1; 
} 

my $report = shift @ARGV;
my $seq2taxid = shift @ARGV;    ## KrakenUniq comes with a file that lists chr<tab>taxID<tab>RefSeq ID + name
my $genbank  = shift @ARGV;     ## Genbank database dump in order to link taxID with representative genome GCF ID
my $bacteria = 0; 

open REPORT,"<",$report or die "$!"; 
open SEQ2TAXID,"<",$seq2taxid or die "$!";
open GENBANK,"<",$genbank or die "$!"; 

my $GB = {};   ## taxid -> rep. genome GCF
my $ST = {};   ## taxid -> either a uniq GCF, an undef, or a "MANY" (this is when we rely on GB) 

while (<GENBANK>) {
	if (m/\trepresentative genome\t(\d+)\t(\d+)\t.*(GCF_.*?)\t/) {
		my $taxid1 = $1; 
		my $taxid2 = $2;
		my $gcf = $3;
		$GB->{$taxid1} = $gcf; 
		$GB->{$taxid2} = $gcf; 
	} elsif (m/\treference genome\t(\d+)\t(\d+)\t.*(GCF_.*?)\t/) {
		my $taxid1 = $1; 
		my $taxid2 = $2;
		my $gcf = $3;
		$GB->{$taxid1} = $gcf; 
		$GB->{$taxid2} = $gcf; 
	} 
} 

while (<SEQ2TAXID>) { 
	m/\t(\d+)\t(GCF.*?) /;
	my $taxid = $1; 
	my $gcf = $2; 
	if (defined $ST->{$taxid} && $ST->{$taxid}->{gcf} ne $gcf) {
		$ST->{$taxid}->{type} = "MANY";
		$ST->{$taxid}->{gcf} = $gcf;  ## keep overwriting - does not matter
	} else { 
		$ST->{$taxid}->{type} = "ONE";
		$ST->{$taxid}->{gcf} = $gcf;
	} 
}

while (<REPORT>) {
	if (m/\t2\tsuperkingdom\t\s+Bacteria/) { 
		$bacteria = 1;
	} elsif ($bacteria == 1) {
		last if (m/\tsuperkingdom\t/);  ## if you've reached the next superkingdom, you're done
		m/\t(\d+)\t([\w\s]+)\t\s+/; 
		my $taxid = $1;
		my $rank = $2;
		if (defined $ST->{$taxid}->{gcf} && $ST->{$taxid}->{type} eq "ONE") { 
			printf "%s\n",$ST->{$taxid}->{gcf};
			printf STDERR "Taxid %s (rank %s): found unique matching sequence ID %s\n",$taxid,$rank,$ST->{$taxid}->{gcf}; 
		} elsif (defined $ST->{$taxid}->{gcf} && $ST->{$taxid}->{type} eq "MANY" && defined $GB->{$taxid}) {
			printf "%s\n",$GB->{$taxid}; 
			printf STDERR "Taxid %s (rank %s): using representative or reference genome sequence ID %s\n",$taxid,$rank,$GB->{$taxid};
		} elsif (defined $ST->{$taxid}->{gcf} && $ST->{$taxid}->{type} eq "MANY" && ! defined $GB->{$taxid}) {
			printf "%s\n",$ST->{$taxid}->{gcf};
			printf STDERR "WARNING: Taxid %s (rank %s) is matched to many GCF IDs, but no representative genome could be found! Using genome %s\n",$taxid,$rank,$ST->{$taxid}->{gcf};
    }
	}
}

close REPORT;
close SEQ2TAXID; 
close GENBANK;

