#!/bin/bash

## v3 of hisat2-based metagenome script 
## - use special script to convert any single cell data into single-end reads with barcodes in name, for umitools; 
## - use KrakenUniq for metagenomic classification and hence have an exact GCF for all used genomes; 
## - use bbduk.sh trimming to trim TSO and other adapters 
## - use slightly more permissive Hisat2 mapping options? (maybe not needed if trimmed well)
## - remove human rRNA reads for speed
## - use igraph to parse the read connectivity graph 
## - use umitools to generate the final output 
TAG=$1

CPUS=4 
SOLODIR=/lustre/scratch126/cellgen/cellgeni/Metagenome/Raquel_metagenome/STARsolo
NCBI=/nfs/cellgeni/NCBI_datasets/bacteria_ref/ncbi_dataset/data/
KUDB=/nfs/cellgeni/Kraken/KrakenUniq_select
SEQ2TAXID=$KUDB/seqid2taxid.map.refseq
GENBANK=$KUDB/assembly_summary_bacteria.txt
GSA=$KUDB/gcf_species_genus.tsv
FCDIR=/nfs/users/nfs_c/cellgeni-su/subread/bin/

#mkdir $TAG
cd $TAG

cp ../*.pl ../*.sh ../parse_read_network.R .
## step 1: make single end fastq file specially formatted for UMI tools
#echo "STEP 1: Making UMI-tools fomatted single-end fastq file .." 
#./make_umi_fastq.pl $SOLODIR/$TAG > se_umi.fastq

## step 1a: need adapter trimming with bbduk - however it's tricky and needs additional evaluation 

## step 2: run Krakenuniq
#echo "STEP 2: Running KrakenUniq .." 
#krakenuniq --threads 16 --db $KUDB --preload-size 120G --output kuniq.output.txt --report-file kuniq.report.txt se_umi.fastq &> kraken_uniq.log 

## step 3: find all RefSeq accessions 
## new script will extract all "precise" matches, and will only keep 1 (representative/reference) genome when many GCFs map to 1 taxid
#echo "STEP 3: Extracting unique bacterial RefSeq accessions .." 
#./kraken_to_GCF.pl kuniq.report.txt $SEQ2TAXID $GENBANK 2> kraken_to_GCF.log | sort | uniq > accessions.list

## remove super-homopolymeric reads - stretches of > 50 of same nucleotide
#echo "STEP 3X: remove super-homopolymeric reads - stretches of > 50 of same nucleotide .." 
#awk '{if (NR%4!=0) {printf "%s\t",$1} else {print $1}}' se_umi.fastq | perl -ne 'print if (! m/A{50}/ && ! m/C{50}/ && ! m/T{50}/ && ! m/G{50}/)' | tr '\t' '\n' > se_umi_filt.fastq

## step 4: now, make combined genome and GFF files for all strains identified by KrakenUniq
## remove the old versions in case they exist
#echo "STEP 4: Generating a combined genome and GFF file for all bacterial species .." 
#if [[ -s combined_bacterial.fna  || -s combined_bacterial.gff ]]
#then
#	rm combined_bacterial.fna combined_bacterial.gff
#fi 

#COUNT=1
#for i in `cat accessions.list`
#do
#  echo "Processing species # $COUNT, RefSeq ID $i .."
#  zcat $NCBI/$i/*genomic.fna.gz >> combined_bacterial.fna & 
#  zcat $NCBI/$i/*genomic.gff.gz | grep -v "^#" >> combined_bacterial.gff &
#  wait 
#  COUNT=$((COUNT+1)) 
#done

## step 5: make hisat2 index
#echo "STEP 5: Making combined Hisat2 index .." 
#hisat2-build combined_bacterial.fna combined_bacterial &> hisat2_build.log

## step 6: align to the produced index. For now we use relaxed soft-clipping penalty; this will disappear with proper adapter trimming
#echo "STEP 6: Mapping original unmapped reads to the combined reference using Hisat2 .." 
#hisat2 -p16 -k 5000 --sp 1,1 --no-spliced-alignment -U se_umi_filt.fastq -x combined_bacterial 2> hisat2_map.log | samtools view -@16 -F4 -O BAM - | samtools sort -@16 -O BAM - > combined_bacterial.bam
#samtools index -@16 combined_bacterial.bam

### step 7: featureCounts it. Multimappers count for 1/N-th of a read. featureCounts compiled in a special way
#echo "STEP 7: Running featureCounts with accurate multimapper counting .." 
#$FCDIR/featureCounts -T 16 -M -O --fraction -t gene -g locus_tag -s 0 -a combined_bacterial.gff -o fcounts.tsv combined_bacterial.bam

### step 8: make a big table of gene biotypes; iterative regex accounts for somewhat complex range of cases that happen in GFF
#echo "STEP 8: Making a gene biotype table .." 
#perl -ne 'if (m/\tgene\t.*.*gene_biotype=(.*?);/) {$type=$1; $lt=""; if (m/;locus_tag=(.*?);/) {$lt=$1} elsif (m/;locus_tag=(.*?)$/) {$lt=$1}; print "$lt\t$type\n"}' combined_bacterial.gff > gene_biotypes.tsv

### step 9: annotate the table and generate a summary per-species 
#echo "STEP 9: Filtering/annotating the featureCounts table and generating summary per RefSeq accession .." 
#./annotate_fc_table.pl fcounts.tsv gene_biotypes.tsv $SEQ2TAXID

## step 10: filter the annotated table based on the number of rRNA and non-rRNA genes; currently require >= 5 reads in each class
#echo "STEP 10: Filtering the accessions, retaining >= 1 rRNA genes, >= 3 non-rRNA genes, and >= 5 reads in each category.." 
#awk -F '\t' '$4>=1 && $5>=5 && $6>=3 && $7 >=5' accession.summary.tsv > filtered.summary.tsv

## step 11: make the network file
#echo "STEP 11: Make the network file using the overlapping reads .." 
#cut -f1 filtered.summary.tsv > filtered.acc.list
#grep -wF -f filtered.acc.list $SEQ2TAXID | cut -f1 > filtered.chr.list
#samtools view -@$CPUS combined_bacterial.bam | grep -wF -f filtered.chr.list | cut -f1,3 | sort | uniq > read_to_chr_filtered_uniq.tsv

#./make_read_network.pl read_to_chr_filtered_uniq.tsv $SEQ2TAXID > nodes_and_edges.tsv

#echo "STEP 12: Parse the network and get the definitive list of strains .." 

#./parse_read_network.R $GSA 0.6

echo "STEP 13: Generate barcode-correction list and correct the barcodes in raw reads .." 
./get_barcode_corrections.sh $TAG $SOLODIR
./correct_barcodes.pl se_umi_filt.fastq barcodes_vs_whitelist.tsv > se_umi_corr.fastq 

echo "STEP 14: Making the new small bacterial reference, and mapping corrected reads to it .."
cut -f1 top_strains.tsv | grep -v RefSeq > top.acc.list
./make_top_reference.sh $NCBI top.acc.list 
hisat2-build top_bacterial.fna top_bacterial &> hisat2_build_top.log 
hisat2 -p$CPUS -k 100 --sp 1,1 --no-spliced-alignment -U se_umi_corr.fastq -x top_bacterial 2> hisat2_top.log | samtools view -@$CPUS -F4 -O BAM - | samtools sort -@$CPUS -O BAM - > top_bacterial.bam
samtools index top_bacterial.bam

echo "STEP 15: Generate a new BAM file with reads assigned to genes using featureCounts .."
## 1) no strand specificity for now; 2) -M -O but without splitting - no point (UMIs will count later) 
$FCDIR/featureCounts -T $CPUS -a top_bacterial.gff -t gene -g locus_tag -M -O -o gene_assigned -R BAM top_bacterial.bam &> top_bacterial.fcounts.log
samtools sort -@$CPUS top_bacterial.bam.featureCounts.bam > top_assigned_sorted.bam
rm top_bacterial.bam.featureCounts.bam
samtools index top_assigned_sorted.bam

echo "STEP 16: Running final quantification using UMI-tools .."
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I top_assigned_sorted.bam -S bac_umi_counts.tsv.gz &> umi_tools_count.log 
