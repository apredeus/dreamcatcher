#!/bin/bash 

CPUS=$1
RAM=$2
SORTCMD="sort -S $((RAM/CPUS))G --parallel=$CPUS"
SORTCMD2="sort -S $((RAM/CPUS/4))G --parallel=$CPUS"

## here we extract unqiue and shared reads, to make networks files with all/rRNA-only/non-rRNA-only reads
## networks are done with shared reads only  
## using the neat efficient parallel GNU sort - matters for huge files 

grep  -F -f filtered.gene.list detected.read_to_gene_w_mismatches.tsv | cut -f1,2 | $SORTCMD | uniq > filtered.read_to_gene_uniq.tsv

$SDIR/annotate_filtered_reads.pl filtered.annotated_fcounts.tsv filtered.read_to_gene_uniq.tsv > filtered.read_type_strain.tsv
rm filtered.read_to_gene_uniq.tsv

cut -f1,4 filtered.read_type_strain.tsv | $SORTCMD | uniq | cut -f1 | uniq -c | awk '{print $2"\t"$1}' > reads_by_strain_uniqness.tsv
awk '$2 == 1' reads_by_strain_uniqness.tsv | cut -f1 | $SORTCMD > filtered.uniq_stain_reads.list
awk '$2 >  1' reads_by_strain_uniqness.tsv | cut -f1 | $SORTCMD > filtered.multi_strain_reads.list

grep -wF -f filtered.uniq_stain_reads.list   filtered.read_type_strain.tsv | awk '$3 == "rRNA"' | cut -f1,4 | $SORTCMD2 | uniq | cut -f2 | $SORTCMD2 | uniq -c | awk '{print $2"\t"$1}' > rRNA.uniq_counts.tsv &
grep -wF -f filtered.uniq_stain_reads.list   filtered.read_type_strain.tsv | awk '$3 != "rRNA"' | cut -f1,4 | $SORTCMD2 | uniq | cut -f2 | $SORTCMD2 | uniq -c | awk '{print $2"\t"$1}' > protein.uniq_counts.tsv &
grep -wF -f filtered.multi_strain_reads.list filtered.read_type_strain.tsv | awk '$3 == "rRNA"' | cut -f1,2 | $SORTCMD2 | uniq > rRNA.shared_reads.tsv &
grep -wF -f filtered.multi_strain_reads.list filtered.read_type_strain.tsv | awk '$3 != "rRNA"' | cut -f1,2 | $SORTCMD2 | uniq > protein.shared_reads.tsv &
wait
