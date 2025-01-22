#!/bin/bash 

NCBI=$1
LIST=$2
TAG=$3

if [[ -e $TAG.fna  || -e $TAG.gff ]]
then
	echo "WARNING: found existing files $TAG.fna/$TAG.gff, deleting and re-generating them .." 
  rm $TAG.fna $TAG.gff
fi  

COUNT=1 
for i in `cat $LIST` 
do 
  echo "Processing species # $COUNT, RefSeq ID $i .."
  zcat $NCBI/$i/*genomic.fna.gz >> $TAG.fna & 
  zcat $NCBI/$i/*genomic.gff.gz | grep -v "^#" >> $TAG.gff &
  wait 
  COUNT=$((COUNT+1)) 
done
