#!/bin/bash 

TAG=$1
CPUS=$2

## this is the master PER SAMPLE stat script 
## three main outputs are: 
## (1) number of reads that maps to detected/filtered/top strains; 
## (2) number of accessions in detected/filtered/top strains; 
## (3) number of protein coding and rRNA genes in detected/filtered/top strains 

## get total, unmapped, unmapped-filtered, and assigned by KrakenUniq as bacterial
RTOT=`perl -ne 'print "$1" if (m/(\d+) total reads, (\d+) reads not mapped/)' 0_preprocessed_reads/preprocess_reads.log`
RUNM=`perl -ne 'print "$2" if (m/(\d+) total reads, (\d+) reads not mapped/)' 0_preprocessed_reads/preprocess_reads.log`
RFLT=`grep "reads; of these"         2_detected_strains/detected.hisat2_map.log | awk '{print $1}'`
RKUQ=`grep "superkingdom.*Bacteria$" 1_krakenuniq/kuniq.report.txt              | cut -f2`

## reads mapped to genomes of detected strains: 
RD1=`grep "aligned exactly 1 time"   2_detected_strains/detected.hisat2_map.log | awk '{print $1}'`
RD2=`grep "aligned >1 times"         2_detected_strains/detected.hisat2_map.log | awk '{print $1}'`
RD3=$((RD1+RD2))

## now, how many of the mapped reads were assigned to a (RefSeq) gene; also, how many were rRNA and non-rRNA
## (note that tRNA/tmRNA genes etc are included into "protein" genes for convenience
RD4=`grep -F -f 2_detected_strains/detected.protein_genes.list 2_detected_strains/detected.read_to_gene_w_mismatches.tsv | cut -f1 | sort | uniq | wc -l`
RD5=`grep -F -f 2_detected_strains/detected.rRNA_genes.list    2_detected_strains/detected.read_to_gene_w_mismatches.tsv | cut -f1 | sort | uniq | wc -l`
RD6=$((RD4+RD5))

## another simple estimate: total detected genes, and rRNA/protein coding breakdown
GD1=`cat 2_detected_strains/detected.protein_genes.list | wc -l`
GD2=`cat 2_detected_strains/detected.rRNA_genes.list    | wc -l`
GD3=$((GD1+GD2))

## now do accessions, out out necessity - if no filtered strains were produced, all other stats have to be 0
ACK=`awk '/superkingdom\t    Bacteria/,/superkingdom\t    Archaea/'  1_krakenuniq/kuniq.report.txt | grep -w assembly | wc -l`
ACD=`cat 2_detected_strains/detected_strains.list | wc -l`
ACA=`cut -f12 detected.annotated_fcounts.tsv | sort | uniq | wc -l` ## detected + at least 1 annotated gene is covered
ACF=`cat 3_filtered_strains/filtered_strains.list | wc -l`
ACT=`cat 5_top_strains/top_strains.list | wc -l`

## now a key fork. if N(filt)>0, then N(top)>0 - it's a must. RF1..RF6 are similar read counts for filtered strains 
## as definted above for detected strains (RD1..RD6);
## similarly, RT1..RT6 are read counts for top strains; GF1..GF3 are gene counts for filtered, and GT1..GT3 for tops. 
RF1=0
RF2=0
RF3=0
RF4=0
RF5=0
RF6=0

RT1=0
RT2=0
RT3=0
RT4=0
RT5=0
RT6=0

GF1=0
GF2=0
GF3=0
GT1=0
GT2=0
GT3=0

## there is a (slight) difference between filtered genes and genes from filtered strains
## this might become important for some artifact hunting
## so we use filtered.gene.list which has filtered genes that also belong to filtered strains
if [[ $ACF != 0 ]]
then
    ## we don't map to filtered, so this is a replacement
    cut -f2 3_filtered_strains/filtered.annotated_fcounts.tsv | sort | uniq > 3_filtered_strains/filtered.chr.list
    RF1=`$CMD samtools view -@$CPUS 2_detected_strains/detected_fcounts.bam | grep -wF -f 3_filtered_strains/filtered.chr.list | grep  -P "\tNH:i:1\t" | cut -f1 | sort | uniq | wc -l`
    RF2=`$CMD samtools view -@$CPUS 2_detected_strains/detected_fcounts.bam | grep -wF -f 3_filtered_strains/filtered.chr.list | grep -vP "\tNH:i:1\t" | cut -f1 | sort | uniq | wc -l`
    RF3=$((RF1+RF2)) 
    RF4=`grep -vP "\trRNA\t" 3_filtered_strains/filtered.read_type_strain.tsv | cut -f1 | sort | uniq | wc -l`
    RF5=`grep  -P "\trRNA\t" 3_filtered_strains/filtered.read_type_strain.tsv | cut -f1 | sort | uniq | wc -l`
    RF6=$((RF4+RF5))

    ## we mapped to top strains, so this is the same as for detected
    RT1=`grep "aligned exactly 1 time" 5_top_strains/top.hisat2_map.log | awk '{print $1}'`
    RT2=`grep "aligned >1 times"       5_top_strains/top.hisat2_map.log | awk '{print $1}'`
    RT3=$((RT1+RT2))
    RT4=`samtools view -@$CPUS 5_top_strains/top_filtered_gene_sorted.bam | grep -F -f 5_top_strains/top.protein_genes.list | cut -f1 | sort | uniq | wc -l`
    RT5=`samtools view -@$CPUS 5_top_strains/top_filtered_gene_sorted.bam | grep -F -f 5_top_strains/top.rRNA_genes.list    | cut -f1 | sort | uniq | wc -l`
    RT6=$((RT4+RT5))

    ## now count genes for filtered and top strains 
    GF1=`cat 3_filtered_strains/filtered.protein_genes.list | wc -l`
    GF2=`cat 3_filtered_strains/filtered.rRNA_genes.list    | wc -l`
    GF3=$((GF1+GF2))
    GT1=`cat 5_top_strains/top.protein_genes.list | wc -l`
    GT2=`cat 5_top_strains/top.rRNA_genes.list    | wc -l`
    GT3=$((GT1+GT2))
fi

printf "Sample\tTot_reads\tHost_unmapped\tFilt_unmapped\tKU_bac\t"
printf "Strains_KU\tStrains_det\tStrains_gene\tStrains_filt\tStrains_top\t"
printf "Genes_det_prot\tGenes_det_rRNA\tGenes_det_all\t"
printf "Genes_filt_prot\tGenes_filt_rRNA\tGenes_filt_all\t"
printf "Genes_top_prot\tGenes_top_rRNA\tGenes_top_all\t"
printf "Reads_det_uniq\tReads_det_mult\tReads_det_map\tReads_det_prot\tReads_det_rRNA\tReads_det_gene\t"
printf "Reads_filt_uniq\tReads_filt_mult\tReads_filt_map\tReads_filt_prot\tReads_filt_rRNA\tReads_filt_gene\t"
printf "Reads_top_uniq\tReads_top_mult\tReads_top_map\tReads_top_prot\tReads_top_rRNA\tReads_top_gene\n"

printf "%s\t%s\t%s\t%s\t%s\t"      $TAG $RTOT $RUNM $RFLT $RKUQ
printf "%s\t%s\t%s\t%s\t%s\t"      $ACK $ACD  $ACA  $ACF  $ACT
printf "%s\t%s\t%s\t"              $GD1 $GD2  $GD3
printf "%s\t%s\t%s\t"              $GF1 $GF2  $GF3
printf "%s\t%s\t%s\t"              $GT1 $GT2  $GT3
printf "%s\t%s\t%s\t%s\t%s\t%s\t"  $RD1 $RD2  $RD3  $RD4  $RD5  $RD6
printf "%s\t%s\t%s\t%s\t%s\t%s\t"  $RF1 $RF2  $RF3  $RF4  $RF5  $RF6
printf "%s\t%s\t%s\t%s\t%s\t%s\n"  $RT1 $RT2  $RT3  $RT4  $RT5  $RT6
