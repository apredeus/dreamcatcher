#!/bin/bash 

LIST=$1

if [[ $LIST == "" ]]
then 
	echo "Usage: ./find_weird_cases.sh <sample_list>" 
	exit 1
fi

for i in `cat $LIST`
do
	N1=`cat $i/filtered.summary.tsv | wc -l`
  if (( $N1 == 0 ))
	then 
		cp $i/filtered.summary.tsv $i/top_strains.tsv
	fi 
	N2=`grep -v Species_taxid $i/top_strains.tsv | awk -F '\t' '{sum+=$12} END {print sum}'`
	if [[ $N2 == "" ]]
	then 
		N2=0
	fi

	if [[ $N1 != $N2 ]]
	then
		echo "Sample $i: number of strains in filtered.summary.tsv ($N1) and in top_strains.tsv ($N2) does not match!" 
	fi 
done 
