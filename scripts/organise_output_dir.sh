#!/bin/bash 

BULK=$1

if [[ $BULK == false ]] 
then
    mkdir -p 0_preprocessed_reads 1_krakenuniq 2_detected_strains 3_filtered_strains 4_read_networks 5_top_strains 6_umitools
    mv host.bam.flagstat barcode_bowtie2.log barcodes_to_correct.fa barcodes_vs_whitelist.tsv whitelist.fa Unmapped_filt.R?.fastq.gz 0_preprocessed_reads/
    mv kuniq* 1_krakenuniq/
    mv detected* gene_biotypes.tsv 2_detected_strains/
    mv filter* blacklisted_genes.list hisat2_human.log human_remap.bam read_* 3_filtered_strains/
    mv rRNA.* protein.* nodes_and_edges.tsv multi_filt_reads.list unique_filt_reads.list reads_by_uniqness.tsv 4_read_networks/
    mv top* parse_read_network.log 5_top_strains/
    mv gene.umi_counts.tsv.gz strain.umi_counts.tsv.gz gene.umi_tools_count.log strain.umi_tools_count.log barcode_to_species.tsv 6_umitools/

    ## keep key files in the root dir
    cp 1_krakenuniq/kuniq.report.txt .
    cp 2_detected_strains/detected.summary.tsv . 
    cp 2_detected_strains/detected.annotated_fcounts.tsv .
    cp 3_filtered_strains/filtered.annotated_fcounts.tsv .
    cp 3_filtered_strains/filtered.summary.tsv . 
    cp 3_filtered_strains/filtered.cluster.tsv . 
    cp 5_top_strains/top.cluster.tsv . 
    cp 5_top_strains/top.annotated_fcounts.tsv . 
    cp 5_top_strains/top.summary.tsv . 
    cp 6_umitools/barcode_to_species.tsv . 
else
    mkdir -p 0_preprocessed_reads 1_krakenuniq 2_detected_strains 3_filtered_strains 4_read_networks 5_top_strains
    mv host.bam.flagstat Unmapped_filt.R?.fastq.gz 0_preprocessed_reads/
    mv kuniq* 1_krakenuniq/
    mv detected* gene_biotypes.tsv 2_detected_strains/
    mv filter* blacklisted_genes.list hisat2_human.log human_remap.bam read_* 3_filtered_strains/
    mv rRNA.* protein.* nodes_and_edges.tsv multi_filt_reads.list unique_filt_reads.list reads_by_uniqness.tsv 4_read_networks/
    mv top* parse_read_network.log 5_top_strains/

    ## keep key files in the root dir
    cp 1_krakenuniq/kuniq.report.txt .
    cp 2_detected_strains/detected.summary.tsv . 
    cp 2_detected_strains/detected.annotated_fcounts.tsv .
    cp 3_filtered_strains/filtered.annotated_fcounts.tsv .
    cp 3_filtered_strains/filtered.summary.tsv . 
    cp 3_filtered_strains/filtered.cluster.tsv . 
    cp 5_top_strains/top.cluster.tsv . 
    cp 5_top_strains/top.annotated_fcounts.tsv . 
    cp 5_top_strains/top.summary.tsv . 
fi
