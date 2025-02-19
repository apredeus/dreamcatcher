#!/usr/bin/env perl 

## v3 of the script: specific KrakenUniq report format (with pseudo-taxIDs assigned to assemblies and sequences)
## is used to defined the list of assemblies that would cover the results most comprehensive coverage 
## 1) there is a unique k-mer cutoff for genus and species of interest; change those to change the number of the output sequences 
## 2) if there is a genus->species->assembly(ies), choose one with most k-mers for each species passing the cutoff; 
## 3) if there is a genus->species, pick the assembly that's best technically (complete > representative > all others); 
## 4) is there is just a genus, pick one species with a good assembly (it's probably 16S and thus does not matter) 

## this script requires a gcf_species_genus file, that is generated separately from the taxDB using parse_taxdb.pl script 

use strict; 
use warnings; 

if (scalar @ARGV != 3) { 
	print "USAGE: ./kuniq_to_GCF.pl <KrakenUniq_report> <gcf_species_genus_table> <genbank_acc_file>\n"; 
	exit 1; 
} 

my $report = shift @ARGV;
my $gcf2taxid = shift @ARGV;    ## 5-column file generated from taxDB: RefSeq ID (== GCF), species taxid, species name, genus taxid, genus name
my $genbank  = shift @ARGV;     ## Genbank database dump to annotate assembly type
my $bacteria = 0; 
my $genus_cutoff = 10; 
my $species_cutoff = 3; 

open REPORT,"<",$report or die "$!"; 
open GCF2TAXID,"<",$gcf2taxid or die "$!";
open GENBANK,"<",$genbank or die "$!"; 

my $GSA = {};  ## genome-species-assembly tree
my $META = {}; ## flat hash of relevant genus, species, and assembly stats/titles

while (<GENBANK>) {
	my @t = split /\t/;
	my $gcf = $t[0];
	$META->{$gcf}->{st} = $t[6]; 
	$META->{$gcf}->{assembly_level} = $t[11]; 
	$META->{$gcf}->{refseq_category} = $t[4];
	my $assembly_score = 0; 
    if ($META->{$gcf}->{assembly_level} eq "Complete Genome") {
        $assembly_score = 10; 
	} elsif ($META->{$gcf}->{assembly_level} eq "Chromosome") {
        $assembly_score = 8; 
	} elsif ($META->{$gcf}->{assembly_level} eq "Scaffold") { 
		$assembly_score = 2; 
	} elsif ($META->{$gcf}->{assembly_level} eq "Contig") { 
		$assembly_score = 1;
	} 
							
	if ($META->{$gcf}->{refseq_category} eq "reference genome") { 
        $assembly_score += 5; 
	} elsif ($META->{$gcf}->{refseq_category} eq "representative genome") { 
        $assembly_score += 3;
	}
    $META->{$gcf}->{assembly_score} = $assembly_score; 
} 

while (<GCF2TAXID>) {
	chomp; 
	my @t = split /\t/; 
	my $gcf = $t[0];
	my $st = $t[1];
	my $gt = $t[3];
	if ($st ne "NONE" && $gt ne "NONE") {
		printf STDERR "kuniq_to_GCF.pl: WARNING: Species taxid according to Genbank (%d) and according to NCBI taxDB (%d) differ!\n",$META->{$gcf}->{st},$st if ($META->{$gcf}->{st} != $st); 
		push @{$GSA->{$gt}->{$st}},$gcf;
		$META->{$gt}->{name} = $t[4]; 
		$META->{$st}->{name} = $t[2]; 
	} 
} 

my ($gt,$st,$genus_spaces,$species_spaces) = (undef) x 4; 

while (<REPORT>) {
	if (m/\t2\tsuperkingdom\t\s+Bacteria/) { 
		$bacteria = 1;
	} elsif ($bacteria == 1 && m/\tsuperkingdom\t/) {
	    last; 
	} elsif ($bacteria == 1) {
		## parse every bacterial line 
	    my @t = split /\t/;
		my $rank = $t[7]; 
		my $kmers = $t[3];
		$t[8] =~ m/^(\s+)/; 
		my $spaces = length $1; 
    
		if ($rank eq "genus") { 
		    $gt = $t[6]; 
			$genus_spaces = $spaces; 
		    if (! defined $META->{$gt}) {
                ## this really shouldn't happen 
			    printf STDERR "kuniq_to_GCF.pl: WARNING: Genus with taxid %d was not defined in NCBI/Genbank files!\n",$gt;
		    } else {
			    $META->{$gt}->{kmers} = $kmers;
            } 
		} 
		if ($rank eq "species" && $spaces > $genus_spaces) { 
			$st = $t[6];
			$species_spaces = $spaces; 
		    if (! defined $META->{$st}) {
                ## this really shouldn't happen 
			    printf STDERR "kuniq_to_GCF.pl: WARNING: Species (%d) / genus(%d) pair was not defined in NCBI/Genbank files!\n",$st,$gt;
            } else {
                $META->{$st}->{kmers} = $kmers;
            } 
        } 
		if ($rank eq "assembly" && $spaces > $species_spaces) { 
			$t[8] =~ m/(GCF_\d+\.\d+)/;
			my $gcf = $1; 
            if (! defined $META->{$gcf}) {
                ## this really shouldn't happen 
		        printf STDERR "kuniq_to_GCF.pl: WARNING: Assembly %s was not defined for species %d / genus %d in NCBI/Genbank files!\n",$gcf,$st,$gt;
		    } else {
		        $META->{$gcf}->{kmers} = $kmers;
            } 
		} 
	}
}

foreach my $gt (keys %{$GSA}) {
    if (defined $META->{$gt}->{kmers} && $META->{$gt}->{kmers} >= $genus_cutoff) { 
	    printf STDERR "Genus %d (%s) DID pass the k-mer threshold - parsing species .. \n",$gt,$META->{$gt}->{name}; 
		## if genus is detected above cutoff, 2 possible scenarios can happen:
		## 1) no species/assemblies are reported, or pass species cutoff; 2) at least 1 species is reported and is above cutoff
		## in case (1) we report a random species with highest quality assembly; in case (2) we report 1 assembly per passing species
		foreach my $st (keys %{$GSA->{$gt}}) {
		    if (defined $META->{$st}->{kmers} && $META->{$st}->{kmers} >= $species_cutoff) {
		        printf STDERR "Species %d (%s) genus %d (%s) DID pass the k-mer threshold - parsing assemblies .. \n",$st,$META->{$st}->{name},$gt,$META->{$gt}->{name};
				## first we check if there are GCFs that map to this species that has some unique k-mers map to it; we choose the GCF with the highest number of k-mers
				my ($best_gcf,$best_assembly) = (undef) x 2;
				my $best_kmers = 0; 
			    foreach my $gcf (@{$GSA->{$gt}->{$st}}) {
					if (defined $META->{$gcf}->{kmers} && $META->{$gcf}->{kmers} >= $best_kmers) {
						$best_gcf = $gcf;
						$best_kmers = $META->{$gcf}->{kmers}; 
					} 
				}
				## second, if no GCF had unique k-mers mapped to it, we choose the best assembly quality
                if ($best_kmers == 0) { 
                    foreach my $gcf (@{$GSA->{$gt}->{$st}}) { 
                        if (defined $META->{$gcf}->{assembly_score} && $META->{$gcf}->{assembly_score} >= $best_assembly) { 
                            $best_gcf = $gcf; 
							$best_assembly = $META->{$gcf}->{assembly_score}; 
					    }
				    } 
					printf STDERR "Best GCF %s for species %d (%s) determined based on assembly quality: level %s, refseq %s, assembly score %d\n",$best_gcf,$st,$META->{$st}->{name},$META->{$best_gcf}->{assembly_level},$META->{$best_gcf}->{refseq_category},$META->{$best_gcf}->{assembly_score};
				} else { 
					printf STDERR "Best GCF %s for species %d (%s) determined based on k-mer count: unique k-mers = %d\n",$best_gcf,$st,$META->{$st}->{name},$best_kmers;
				}
				print "$best_gcf\n"; 
			} else { 
		        printf STDERR "Species %d (%s) genus %d (%s) did not pass the k-mer threshold\n",$st,$META->{$st}->{name},$gt,$META->{$gt}->{name}; 
			}
		} 
	} else { 
        printf STDERR "Genus %d (%s) did not pass the k-mer threshold\n",$gt,$META->{$gt}->{name};
	} 
} 

close REPORT;
close GCF2TAXID; 
close GENBANK;
