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

CPUS=16 
SOLODIR=/lustre/scratch126/cellgen/cellgeni/Metagenome/Raquel_metagenome/STARsolo
NCBI=/nfs/cellgeni/NCBI_datasets/bacteria_ref/ncbi_dataset/data/
KUDB=/nfs/cellgeni/Kraken2/Kraken_uniq
GENBANK=/nfs/cellgeni/Kraken2/assembly_summary_genbank.txt
FCDIR=/nfs/users/nfs_c/cellgeni-su/subread/bin/

#mkdir $TAG
cd $TAG

cp ../make_umi_fastq.pl ../kraken_to_GCF.pl ../annotate_fc_table.pl .
## step 1: make single end fastq file specially formatted for UMI tools
#echo "STEP 1: Making UMI-tools fomatted single-end fastq file .." 
#./make_umi_fastq.pl $SOLODIR/$TAG > $TAG.se_umi.fastq

## step 1a: need adapter trimming with bbduk - however it's tricky and needs additional evaluation 

## step 2: run Krakenuniq
#echo "STEP 2: Running KrakenUniq .." 
#krakenuniq --threads 16 --db $KUDB --preload-size 120G --output $TAG.kuniq.output.txt --report-file $TAG.kuniq.report.txt $TAG.se_umi.fastq &> kraken_uniq.log 

## step 3: find all RefSeq accessions 
## new script will extract all "precise" matches, and will only keep 1 (representative/reference) genome when many GCFs map to 1 taxid
#echo "STEP 3: Extracting unique bacterial RefSeq accessions .." 
#./kraken_to_GCF.pl $TAG.kuniq.report.txt $KUDB/seqid2taxid.map $GENBANK 2> kraken_to_GCF.log | sort | uniq > accessions.list

## step 4: now, make combined genome and GFF files for all strains identified by KrakenUniq
## remove the old versions in case they exist
#echo "STEP 4: Generating a combined genome and GFF file for all bacterial species .." 
#if [[ -s $TAG.combined.fna  || -s $TAG.combined.gff ]]
#then
#	rm $TAG.combined.fna $TAG.combined.gff
#fi 

#COUNT=1
#for i in `cat accessions.list`
#do
#  echo "Processing species # $COUNT, RefSeq ID $i .."
#  zcat $NCBI/$i/*genomic.fna.gz >> $TAG.combined.fna & 
#  zcat $NCBI/$i/*genomic.gff.gz | grep -v "^#" >> $TAG.combined.gff &
#  wait 
#  COUNT=$((COUNT+1)) 
#done

## step 5: make hisat2 index
#echo "STEP 5: Making combined Hisat2 index .." 
#hisat2-build $TAG.combined.fna $TAG.combined &> hisat2_build.log

## step 6: align to the produced index. For now we use relaxed soft-clipping penalty; this will disappear with proper adapter trimming
echo "STEP 6: Mapping original unmapped reads to the combined reference using Hisat2 .." 
hisat2 -p16 -k 10000 --sp 1,1 --no-spliced-alignment -U $TAG.se_umi.fastq -x $TAG.combined 2> hisat2_map.log | samtools view -@16 -F4 -O BAM - | samtools sort -@16 -O BAM - > $TAG.combined.bam
samtools index -@16 $TAG.combined.bam

## step 7: featureCounts it. Multimappers count for 1/N-th of a read. featureCounts compiled in a special way
echo "STEP 7: Running featureCounts with accurate multimapper counting .." 
$FCDIR/featureCounts -T 16 -M -O --fraction -t gene -g locus_tag -s 0 -a $TAG.combined.gff -o $TAG.fcounts.tsv $TAG.combined.bam

## step 8: make a big table of gene biotypes; iterative regex accounts for somewhat complex range of cases that happen in GFF
echo "STEP 8: Making a gene biotype table .." 
perl -ne 'if (m/\tgene\t.*.*gene_biotype=(.*?);/) {$type=$1; $lt=""; if (m/;locus_tag=(.*?);/) {$lt=$1} elsif (m/;locus_tag=(.*?)$/) {$lt=$1}; print "$lt\t$type\n"}' $TAG.combined.gff > gene_biotypes.tsv

## step 9: annotate the table and generate a summary per-species 
echo "STEP 9: Filtering/annotating the featureCounts table and generating summary per RefSeq accession .." 
./annotate_fc_table.pl $TAG.fcounts.tsv gene_biotypes.tsv $KUDB/seqid2taxid.map
