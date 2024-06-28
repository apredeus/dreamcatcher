#!/bin/bash 

SMP=$1
SOLODIR=$2
if [[ $SMP == "" || $SOLODIR == "" ]]
then
	>&2 echo "Usage: ./get_barcode_corrections.sh <sample_id> <starsolo_dir>" 
	exit 1
fi 

CPUS=4

echo "Generating a FASTA file from unique cellular barcodes that need to be corrected against the whitelist.." 
awk 'NR%4==1' se_umi_filt.fastq | perl -ne '@t=split/_/; print "$t[-2]\n"' | sort | uniq | awk '{print ">"$1; print $1}' > se_umi_filt.barcodes.fa

WL=`grep "soloCBwhitelist" $SOLODIR/$SMP/Log.out  | grep RE-DEFINED | awk '{print $2}'`
echo "Found whitelist $WL; making a fasta file .." 
cat $WL | awk '{print ">"$1; print $1}' > whitelist.fa
bowtie2-build whitelist.fa whitelist &> /dev/null 

echo "Mapping barcodes to the whitelist .." 
bowtie2 -af -p$CPUS --no-hd --very-sensitive -U se_umi_filt.barcodes.fa -x whitelist > barcodes_vs_whitelist.sam 2> barcode_bowtie2.log
perl -ne '@t=split/\t/; print "$t[0]\t$t[2]\t$1\n" if ($t[5] eq "16M" && ($t[1] == 0 || $t[1] == 256) && m/\tNM:i:([0,1])\t/)' barcodes_vs_whitelist.sam > barcodes_vs_whitelist.tsv
