#!/bin/bash 

if [[ ! -s top.cluster.tsv ]] 
then
  echo "ERROR: top.cluster.tsv does not exist or is empty, no reason to continue this. I quit!"
  exit 1
fi


## generate a list of all prefixes
KK=`zcat bac_umi_counts.tsv.gz | grep -v gene | cut -f1 | tr ',' '\n' | perl -ne 's/\d+$//; print' | sort | uniq`

## some prefixes are bad (on purpose!) - this should fix it

for i in $KK
do 
  for j in $KK
  do 
    if [[ $i != $j && `echo $j | grep "^$i"` != "" ]]
	then
	  ## remove a longer prefix if one is contained inside another
	  KK=`echo $KK | sed "s/$j//"`
	  >&2 echo "WARNING: gene prefixes $i and $j collide (and probably match the same assembly!). Removing $j.." 
	fi
  done
done

for i in $KK
do
  GCF=`grep "^$i" filtered.annotated_fcounts.tsv | cut -f12 | uniq`
  SPECIES=`grep $GCF top.cluster.tsv | cut -f11`
  echo -e "$i\t$SPECIES"
done 
