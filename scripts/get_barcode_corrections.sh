#!/bin/bash 

## this is ONLY necessary for STARsolo Unmapped reads fastq output - BAMs contain corrected barcodes already
## (yes, even for the unmapped reads)

FQDIR=$1
TAG=$2
UNCORR=Unmapped_filt_uncorr.R1.fastq

if [[ $TAG == "" || $FQDIR == "" ]]
then
	>&2 echo "Usage: ./get_barcode_corrections.sh <starsolo_dir> <sample_id>" 
	exit 1
fi 

echo -e "get_barcode_corrections.sh: generating a FASTA file from unique cellular barcodes that need to be corrected against the whitelist.." 
awk 'NR%4==1' $UNCORR | perl -ne '@t=split/_/; print "$t[-2]\n"' | sort | uniq | awk '{print ">"$1; print $1}' > barcodes_to_correct.fa

WL=`grep "soloCBwhitelist" $FQDIR/$TAG/Log.out  | grep RE-DEFINED | awk '{print $2}'`
echo "get_barcode_corrections.sh: found whitelist $WL; making a fasta file .." 
cat $WL | awk '{print ">"$1; print $1}' > whitelist.fa
$CMD bowtie2-build whitelist.fa whitelist &> /dev/null 

echo -e "get_barcode_corrections.sh: mapping barcodes to the whitelist .." 
$CMD bowtie2 -af -p$CPUS --no-hd --very-sensitive -U barcodes_to_correct.fa -x whitelist > barcodes_vs_whitelist.sam 2> barcode_bowtie2.log
perl -ne '@t=split/\t/; print "$t[0]\t$t[2]\t$1\n" if ($t[5] eq "16M" && ($t[1] == 0 || $t[1] == 256) && m/\tNM:i:([0,1])\t/)' barcodes_vs_whitelist.sam > barcodes_vs_whitelist.tsv
