#!/bin/bash 

CPUS=$1

## here we extract unqiue and shared reads, to make networks files with all/rRNA-only/non-rRNA-only reads
## networks are done with shared reads only  
grep  -F -f filtered.gene.list read_to_gene_detected_w_mismatches.tsv | cut -f1,2 | sort -S 50% --parallel=$CPUS  | uniq > read_to_gene_filtered_uniq.tsv
grep  -F -f filtered.gene.list detected.annotated_fcounts.tsv > filtered.annotated_fcounts.tsv

$SDIR/annotate_filtered_reads.pl filtered.annotated_fcounts.tsv read_to_gene_filtered_uniq.tsv > read_filtered_type_strain.tsv
cut -f1,4 read_filtered_type_strain.tsv | sort -S 50% --parallel=$CPUS | uniq | cut -f1 | uniq -c | awk '{print $2"\t"$1}' > reads_by_uniqness.tsv
awk '$2 == 1' reads_by_uniqness.tsv | cut -f1 | sort > unique_filt_reads.list &
awk '$2 >  1' reads_by_uniqness.tsv | cut -f1 | sort > multi_filt_reads.list  &
wait

grep -wF -f unique_filt_reads.list read_filtered_type_strain.tsv | awk '$3 == "rRNA"' | cut -f1,4 | sort | uniq | cut -f2 | sort | uniq -c | awk '{print $2"\t"$1}' > rRNA.uniq_counts.tsv &
grep -wF -f unique_filt_reads.list read_filtered_type_strain.tsv | awk '$3 != "rRNA"' | cut -f1,4 | sort | uniq | cut -f2 | sort | uniq -c | awk '{print $2"\t"$1}' > protein.uniq_counts.tsv &
grep -wF -f multi_filt_reads.list  read_filtered_type_strain.tsv | awk '$3 == "rRNA"' | cut -f1,2 | sort | uniq > rRNA.shared_reads.tsv &
grep -wF -f multi_filt_reads.list  read_filtered_type_strain.tsv | awk '$3 != "rRNA"' | cut -f1,2 | sort | uniq > protein.shared_reads.tsv &
wait
