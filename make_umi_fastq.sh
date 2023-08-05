#!/bin/bash 

DIR=$1 ## could be STARsolo or Cell Ranger dir 
DIR=`readlink -f $DIR`
TAG=`basename $DIR` 

if [[ -d $DIR/output && -s $DIR/Log.final.out && -s $DIR/Unmapped.out.mate1.gz && -s $DIR/Unmapped.out.mate2.gz ]]
then
  >&2 echo "Sample was determined to be STARsolo output! Making a single-end, UMI-tools-formatted fastq file.."
  R1=$DIR/Unmapped.out.mate1.gz
  R2=$DIR/Unmapped.out.mate2.gz
  paste <(zcat $R1) <(zcat $R2)| awk '{if (NR%4!=0) {printf "%s\t",$0} else {print}}' | perl -ne '@t=split/\t/; $rname=$t[0]; $rname=~s/ .*//; $len=length($t[3]); $bclen=($len > 24) ? 16 : 14; $umilen = ($len > 28) ? 10 : $len-$bclen; $bc = substr($t[3],0,$bclen); $umi=substr($t[3],$bclen,$umilen); printf "%s_%s_%s\n%s\n+\n%s\n",$rname,$bc,$umi,$t[2],$t[6]' > $TAG.se_umi.fastq
elif [[ -d $DIR/outs && ( -s $DIR/possorted_genome_bam.bam || -s $DIR/gex_possorted_bam.bam ) ]]
then 
  >&2 echo "Sample was determined to be Cell Ranger output!" 
else 
  >&2 echo "ERROR: the sample could not be matched to neither STARsolo nor Cell Ranger format (or some necessary files are missing)! Exiting.." 
  exit 1
fi 
