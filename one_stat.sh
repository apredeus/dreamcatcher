#!/bin/bash

TAG=$1

## read stats - starting from a total number
## some stats from STARsolo
TOT=`grep "Number of input reads" ../STARsolo/$TAG/Log.final.out | awk '{print $6}'`
UNM=`grep "Number of reads unmapped" ../STARsolo/$TAG/Log.final.out | awk -F "|" '{print $2}' | awk '{sum+=$1} END {print sum}'`

## now how many were filtered (after homopolymer removal), and how many were mapped/classified
FLT=`grep "reads; of these" $TAG/hisat2_map.log | awk '{print $1}'`
KUR=`grep "superkingdom.*Bacteria$" $TAG/kuniq.report.txt | cut -f2`
HSU=`grep "aligned exactly 1 time" $TAG/hisat2_map.log | awk '{print $1}'`
HSM=`grep "aligned >1 times" $TAG/hisat2_map.log | awk '{print $1}'`
HSA=$((HSU+HSM))

## now, how many of the mapped reads were assigned to a (RefSeq) gene; also, how many were rRNA and non-rRNA
grep -P "\trRNA\t" $TAG/annotated.fcounts.tsv | cut -f1 | sort | uniq > $TAG/rRNA_genes.list
grep -P "\tprotein_coding\t" $TAG/annotated.fcounts.tsv | cut -f1 | sort | uniq > $TAG/protein_genes.list
RNA=`grep -F -f $TAG/rRNA_genes.list $TAG/combined_bacterial_fcount_assignment.tsv | cut -f1 | sort | uniq | wc -l`
PRT=`grep -F -f $TAG/protein_genes.list $TAG/combined_bacterial_fcount_assignment.tsv | cut -f1 | sort | uniq | wc -l`
ASS=$((RNA+PRT))

## another simple estimates: total detected genes, and rRNA/protein coding breakdown 
NRNA=`cat $TAG/rRNA_genes.list | wc -l`
NPRT=`cat $TAG/protein_genes.list | wc -l`
NGEN=$((NRNA+NPRT))

## now do accessions, out out necessity - if no filtered strains were produced, all other stats have to be 0
ACK=`awk '/superkingdom\t    Bacteria/,/superkingdom\t    Archaea/'  $TAG/kuniq.report.txt | grep -w assembly | wc -l`
ACC=`cat $TAG/accessions.list | wc -l`
ACA=`cut -f9 $TAG/annotated.fcounts.tsv | sort | uniq | wc -l` ## at least 1 gene detected
ACF=`cat $TAG/filtered.acc.list | wc -l`
ACT=`cat $TAG/top.acc.list | wc -l`

## now a key fork. if N(filt)>0, then N(top)>0 - it's a must. 
FRNA=0
FPRT=0
FASS=0
TRNA=0
TPRT=0
TASS=0
THSA=0
THSU=0
THSM=0

if [[ $ACF != 0 ]]
then
  THSU=`grep "aligned exactly 1 time" $TAG/hisat2_top.log | awk '{print $1}'`
  THSM=`grep "aligned >1 times" $TAG/hisat2_top.log | awk '{print $1}'`
  THSA=$((THSU+THSM))
  grep -P "\trRNA\t" $TAG/annotated.fcounts.tsv | grep -F -f $TAG/filtered.acc.list | cut -f1 | sort | uniq > $TAG/filtered_rRNA_genes.list
  grep -P "\tprotein_coding\t" $TAG/annotated.fcounts.tsv | grep -F -f $TAG/filtered.acc.list | cut -f1 | sort | uniq > $TAG/filtered_protein_genes.list
	FRNA=`grep -F -f $TAG/filtered_rRNA_genes.list $TAG/combined_bacterial_fcount_assignment.tsv | cut -f1 | sort | uniq | wc -l`
	FPRT=`grep -F -f $TAG/filtered_protein_genes.list $TAG/combined_bacterial_fcount_assignment.tsv | cut -f1 | sort | uniq | wc -l`
	FASS=$((FRNA+FPRT))
  TRNA=`samtools view $TAG/top_assigned_sorted.bam | grep -F -f $TAG/rRNA_genes.list | cut -f1 | sort | uniq | wc -l`
  TPRT=`samtools view $TAG/top_assigned_sorted.bam | grep -F -f $TAG/protein_genes.list | cut -f1 | sort | uniq | wc -l`
	TASS=$((TRNA+TPRT))
fi

echo -e "$TAG\t$TOT\t$UNM\t$FLT\t$KUR\t$HSA\t$HSU\t$HSM\t$ASS\t$PRT\t$RNA\t$NGEN\t$NPRT\t$NRNA\t$ACK\t$ACC\t$ACA\t$ACF\t$ACT\t$THSA\t$THSU\t$THSM\t$FASS\t$FPRT\t$FRNA\t$TASS\t$TPRT\t$TRNA" > $TAG.stats.tsv
