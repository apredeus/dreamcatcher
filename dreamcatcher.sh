#!/bin/bash

## v3.3 of hisat2-based metagenome script 
## - use special script to convert any single cell data into single-end reads with barcodes in name, for umitools; 
## - use KrakenUniq for metagenomic classification and hence have an exact GCF for all used genomes; 
## - use slightly more permissive than default Hisat2 mapping options
## - use igraph to parse the read connectivity graph 
## - use umitools to generate the final output 


##########################################################################
######################   conda activate dreamcatcher   ###################
########## (until we overhaul the whole thing to Singularity) ############
##########################################################################

SOLODIR=$1
TAG=$2

if [[ $SOLODIR == "" || $TAG == "" ]]
then 
	>&2 echo "Usage: ./dreamcatcher.sh <starsolo/cellranger dir> <sample_id>"
	exit 1
fi

CPUS=16
NCBI=/nfs/cellgeni/NCBI_datasets/bacteria_ref/ncbi_dataset/data/
KUDB=/nfs/cellgeni/Kraken/KrakenUniq_select
SEQ2TAXID=$KUDB/seqid2taxid.map.refseq
GENBANK=$KUDB/assembly_summary_bacteria.txt
GSA=$KUDB/gcf_species_genus.tsv
FCDIR=/nfs/users/nfs_c/cellgeni-su/subread/bin/

mkdir $TAG
cd $TAG

cp ../*.pl ../*.sh ../parse_read_network.R .
## step 1: make single end fastq file specially formatted for UMI tools
## remove homopolymers and correct the barcodes too 
echo "STEP 1: Making UMI-tools fomatted single-end fastq file .." 
./make_umi_fastq.pl $SOLODIR/$TAG se_umi_filt.fastq
./get_barcode_corrections.sh $TAG $SOLODIR
./correct_barcodes.pl se_umi_filt.fastq barcodes_vs_whitelist.tsv se_umi_corr.fastq 
pigz se_umi_corr.fastq
rm se_umi_filt.fastq barcodes_vs_whitelist.sam whitelist*.bt2

# step 2: run Krakenuniq
echo "STEP 2: Running KrakenUniq .." 
krakenuniq --threads $CPUS --db $KUDB --preload-size 120G --output kuniq.output.txt --report-file kuniq.report.txt se_umi_corr.fastq.gz &> /dev/null
pigz kuniq.output.txt

if (( `cat kuniq.report.txt | wc -l` < 10 ))
then 
	>&2 echo "ERROR: STEP 2 (KrakenUniq) has produced an empty report file! Exiting .." 
	exit 1
fi
## step 3: find all RefSeq accessions 
## new script will extract all "precise" matches, and will only keep 1 (representative/reference) genome when many GCFs map to 1 taxid
echo "STEP 3: Extracting unique bacterial RefSeq accessions .."
./kuniq_to_GCF.pl kuniq.report.txt $GSA $GENBANK > accessions.list 2> kuniq_to_GCF.log

## step 4: now, make combined genome and GFF files for all strains identified by KrakenUniq
## remove the old versions in case they exist
echo "STEP 4: Generating a combined genome and GFF file for all bacterial species .." 
./make_combined_reference.sh $NCBI accessions.list combined_bacterial &> mkref_combined.log  

## step 5: make hisat2 index
echo "STEP 5: Making combined Hisat2 index .." 
hisat2-build combined_bacterial.fna combined_bacterial &> hisat2_build.log

## step 6: align to the produced index. For now we use relaxed soft-clipping penalty; this will disappear with proper adapter trimming
echo "STEP 6: Mapping original unmapped reads to the combined reference using Hisat2 .." 
hisat2 -p$CPUS -k 5000 --sp 1,1 --no-spliced-alignment -U se_umi_corr.fastq.gz -x combined_bacterial 2> hisat2_map.log | samtools view -@$CPUS -F4 -O BAM - | samtools sort -@$CPUS --verbosity 0 -O BAM - > combined_bacterial.bam

### step 7: featureCounts it. Multimappers count for 1/N-th of a read. featureCounts compiled in a special way
echo "STEP 7: Running featureCounts with accurate multimapper counting .." 
$FCDIR/featureCounts -T $CPUS -M -O --fraction -t gene -g locus_tag -s 0 -a combined_bacterial.gff -o fcounts.tsv -R CORE combined_bacterial.bam &> combined_bacterial.fcounts.log  
mv combined_bacterial.bam.featureCounts combined_bacterial_fcount_assignment.tsv

### step 8: make a big table of gene biotypes; iterative regex accounts for somewhat complex range of cases that happen in GFF
echo "STEP 8: Making a gene biotype table .." 
perl -ne 'if (m/\tgene\t.*.*gene_biotype=(.*?);/) {$type=$1; $lt=""; if (m/;locus_tag=(.*?);/) {$lt=$1} elsif (m/;locus_tag=(.*?)$/) {$lt=$1}; print "$lt\t$type\n"}' combined_bacterial.gff > gene_biotypes.tsv

### step 9: annotate the table and generate a summary per-species 
echo "STEP 9: Filtering/annotating the featureCounts table and generating summary per RefSeq accession .." 
./annotate_fc_table.pl fcounts.tsv gene_biotypes.tsv $SEQ2TAXID

## step 10: filter the annotated table based on the number of rRNA and non-rRNA genes; currently require >= 5 reads in each class
echo "STEP 10: Filtering the accessions: keeping strains with 5 or more genes detected, and with 5 or more reads.." 
awk -F '\t' '$4+$6 >= 5 && $5+$7 >= 5' accession.summary.tsv > filtered.summary.tsv

if [[ ! -s filtered.summary.tsv ]]
then
	>&2 echo "ERROR: no filtered strains were selected at STEP 10! Exiting .." 
	exit 1
fi

## step 11: make the network file
echo "STEP 11: Make the network file using the overlapping reads .." 
cut -f1 filtered.summary.tsv > filtered.acc.list
grep -wF -f filtered.acc.list $SEQ2TAXID | cut -f1 > filtered.chr.list
grep -wF -f filtered.acc.list annotated.fcounts.tsv > filtered.fcounts.tsv
cut -f1 filtered.fcounts.tsv | sort | uniq > filtered.gene.list
grep  -F -f filtered.gene.list combined_bacterial_fcount_assignment.tsv | cut -f1,4 | sort | uniq > read_to_gene_filtered_uniq.tsv
#samtools view -@$CPUS combined_bacterial.bam | grep -wF -f filtered.chr.list | cut -f1,3 | sort | uniq > read_to_chr_filtered_uniq.tsv

./make_read_network2.pl read_to_gene_filtered_uniq.tsv filtered.fcounts.tsv > nodes_and_edges.tsv

echo "STEP 12: Parse the network and get the definitive list of strains .." 
./parse_read_network.R $GSA 0.6 &> parse_read_network.log

echo "STEP 13: Making the new small bacterial reference, and mapping corrected reads to it .."
cut -f1 top_strains.tsv | grep -v RefSeq > top.acc.list
./make_combined_reference.sh $NCBI top.acc.list top_bacterial &> mkref_top.log 
hisat2-build top_bacterial.fna top_bacterial &> hisat2_build_top.log 
hisat2 -p$CPUS -k 100 --sp 1,1 --no-spliced-alignment -U se_umi_corr.fastq.gz -x top_bacterial 2> hisat2_top.log | samtools view -@$CPUS -F4 -O BAM - | samtools sort -@$CPUS --verbosity 0 -O BAM - > top_bacterial.bam

echo "STEP 14: Generate a new BAM file with reads assigned to genes using featureCounts .."
## 1) no strand specificity for now; 2) -M -O but without splitting - no point (UMIs will count later) 
$FCDIR/featureCounts -T $CPUS -a top_bacterial.gff -t gene -g locus_tag -M -O -o gene_assigned -R BAM top_bacterial.bam &> top_bacterial.fcounts.log
samtools sort -@$CPUS top_bacterial.bam.featureCounts.bam > top_assigned_sorted.bam
rm top_bacterial.bam.featureCounts.bam
samtools index top_assigned_sorted.bam

echo "STEP 15: Running final quantification using UMI-tools .."
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I top_assigned_sorted.bam -S bac_umi_counts.tsv.gz &> umi_tools_count.log 

echo "STEP 16: Generate barcode-to-species table from per-gene UMI-tools .." 
./make_prefix_table.sh > prefix_to_species.tsv
./barcode_to_species.pl bac_umi_counts.tsv.gz prefix_to_species.tsv > barcode_to_species.tsv

echo "ALL DREAMCATCHER PROCESSING IS DONE! ENJOY!!!" 
