#!/bin/bash 

NCBI=/nfs/cellgeni/NCBI_datasets/bacteria_ref/ncbi_dataset/data
LIST=GCF_for_chrom.list
TX=/lustre/scratch126/cellgen/cellgeni/Metagenome/Raquel_metagenome/align3/taxid_acc_name.txt

for i in `cat $LIST`
do
	CHR=`zcat $NCBI/$i/*genomic.fna.gz | grep -F ">" | awk '{print $1}' | sed "s/>//"`
  for j in $CHR
	do
		grep -F $i $TX | awk -v v=$j '{print v"\t"$0}'
	done
done
