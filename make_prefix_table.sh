#!/bin/bash 

if [[ ! -s top_strains.tsv ]] 
then
  echo "ERROR: top_strains.tsv is empty, no reason to continue this. I quit!"
  exit 1
fi


## generate a list of all prefixes
KK=`zcat bac_umi_counts.tsv.gz | grep -v gene | cut -f1 | tr ',' '\n' | perl -ne 's/\d+$//; print' | sort | uniq`

for i in $KK
do
  GCF=`grep "^$i" annotated.fcounts.tsv | cut -f12 | uniq`
  SPECIES=`grep $GCF top_strains.tsv | cut -f3`
  echo -e "$i\t$SPECIES"
done 
