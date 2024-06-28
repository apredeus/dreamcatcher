#!/bin/bash 

## calculate overlaps between the KrakenUniq and mapping/read counting 

TAG=$1

mkdir $TAG/mapping_stats

## get KUniq reads that are assigned to *anything* remotely bacterial 
awk '/superkingdom\t    Bacteria/,/superkingdom\t    Archaea/'  $TAG/kuniq.report.txt | grep -vw Archaea | cut -f 7 | awk '{print "\t"$1"\t"}' > $TAG/mapping_stats/bac_taxid.txt
zcat $TAG/kuniq.output.txt.gz | grep -f $TAG/mapping_stats/bac_taxid.txt > $TAG/mapping_stats/bac_reads.txt

## separate 10-digit (sequence) taxids and other stuff
awk 'length($1) == 10' bac_taxid.txt > seq_taxid.txt
awk 'length($1) != 10' bac_taxid.txt > noseq_taxid.txt

## make read name files
cut -f1 $TAG/combined_bacterial_fcount_assignment.tsv | sort | uniq > $TAG/mapping_stats/mapped_read_names.list
cut -f2 $TAG/mapping_stats/bac_reads.txt | sort | uniq > $TAG/mapping_stats/bac_read_names.list 

## N1/N2: bac KU, not mapped, not assigned, 10 digit ID/regular NCBI
N1=`grep -vw -f mapped_read_names.list bac_reads.txt | grep -c -f seq_taxid.txt`
N2=`grep -vw -f mapped_read_names.list bac_reads.txt | grep -c -f noseq_taxid.txt`

## N3/N4: bac KU, mapped, not assigned, 10 digit ID/regular NCBI
N3=`grep -w -f mapped_read_names.list bac_reads.txt
