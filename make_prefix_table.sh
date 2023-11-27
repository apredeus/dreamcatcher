#!/bin/bash 

if [[ ! -s top_strains.tsv ]] 
then
	echo "ERROR: top_strains.tsv is empty, no reason to continue this. I quit!"
	exit 1
fi


## generate a list of all prefixes
KK=`zcat bac_umi_counts.tsv.gz | grep -v gene | cut -f1 | perl -ne 's/\d+$//g; print' | sort | uniq | grep -v ","`

for i in $KK
do
	GCF=`grep -F $i annotated.fcounts.tsv | cut -f9 | uniq`
  SPECIES=`grep $GCF top_strains.tsv | cut -f3`
	echo -e "$i\t$SPECIES"
done 
