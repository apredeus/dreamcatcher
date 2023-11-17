#!/bin/bash 

NCBI=$1
LIST=$2

if [[ -e top_bacterial.fna  || -e top_bacterial.gff ]]
then 
  rm top_bacterial.fna top_bacterial.gff
fi  

COUNT=1 
for i in `cat $LIST` 
do 
  echo "Processing species # $COUNT, RefSeq ID $i .."
  zcat $NCBI/$i/*genomic.fna.gz >> top_bacterial.fna & 
  zcat $NCBI/$i/*genomic.gff.gz | grep -v "^#" >> top_bacterial.gff &
  wait 
  COUNT=$((COUNT+1)) 
done
