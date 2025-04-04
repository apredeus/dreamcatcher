#!/bin/bash 

FQDIR=$1
TAG=$2
BULK=$3
BAM=$4

if [[ $# != 4 ]]
then
    >&2 echo "Usage: preprocess_unmapped_reads.sh <fastqs_dir> <sample_id> <bulk_or_not> <bam_or_not>"
    exit 1
fi

## this is a key preprocessing script which should take care of all supported input data cases
## and bring it to a common denominator - either one file (Unmapped_filt.R1.fastq.gz)
## or a pair Unmapped_filt.R1.fastq.gz/Unmapped_filt.R2.fastq.gz
## all 10x is converted to single-end, and BC+UMI sequences are added to the read name separated by _ (UMI-tools format) 
## singularity image $CMD and script dir $SDIR are exported to the $ENV in the master script 

if [[ $BAM == true && $BULK == false ]]
then
    BAMFILE=`find $FQDIR/$TAG/* | grep -iv atac | grep -v vdj | grep -v multi | grep "\.bam$"`
    $CMD samtools flagstat -@4 $BAMFILE > host.bam.flagstat
elif [[ $BAM == true && $BULK == true ]]
then
    BAMFILE=`find $FQDIR/$TAG/* | grep -iv transcriptome | grep "\.bam$"`
    $CMD samtools flagstat -@4 $BAMFILE > host.bam.flagstat
fi

if [[ $BULK == false ]]
then 
    ## no PE reads for single-cell; in case of PE 10x, we only use the "big" read (R2), take BC+UMI from R1, and drop the rest of the read
    if [[ $BAM == false ]] 
    then
        ## case 1: STARsolo unmapped reads. Barocodes need to be corrected. 
        echo -e "\tpreprocess_unmapped_reads.sh: preparing unmapped reads in single cell mode, STARsolo raw unmapped reads are used." 
        $SDIR/make_umi_fastq.pl $FQDIR/$TAG Unmapped_filt_uncorr.R1.fastq ## this script filters homopolymers too 
        $SDIR/get_barcode_corrections.sh $TAG $FQDIR
        $SDIR/correct_barcodes.pl Unmapped_filt_uncorr.R1.fastq barcodes_vs_whitelist.tsv Unmapped_filt.R1.fastq
        $CMD pigz Unmapped_filt.R1.fastq
        rm Unmapped_filt_uncorr.R1.fastq barcodes_vs_whitelist.sam whitelist*.bt2
    else
        ## case 2: STARsolo/Cell Ranger/Space Ranger 10x BAM file
        echo -e "\tpreprocess_unmapped_reads.sh: preparing unmapped reads in single cell mode, 10x BAM file is used." 
        $SDIR/make_umi_fastq.pl $FQDIR/$TAG Unmapped_filt.R1.fastq ## this script filters homopolymers too
        $CMD pigz Unmapped_filt.R1.fastq
    fi
else
    ## for bulk, reads can come in BAM or STAR fastq, and can be SE/PE
    echo -e "\tpreprocess_unmapped_reads.sh: preparing unmapped reads in bulk mode.." 
    $SDIR/make_bulk_fastq.pl $FQDIR/$TAG  ## filter homopolymers and remove reads shorter than 40 bp
    if [[ -s Unmapped_filt.R1.fastq && -s Unmapped_filt.R2.fastq ]]
    then 
        ## cases 3,4: paired-end reads (either from BAM or from unmapped fastq) 
        $CMD pigz Unmapped_filt.R1.fastq
        $CMD pigz Unmapped_filt.R2.fastq
        if [[ -s Unmapped_unfilt.R1.fastq ]] 
        then 
          rm Unmapped_unfilt.*.fastq
        fi
    else
        ## cases 5,6: single-end reads (either from BAM or from unmapped fastq)
        $CMD pigz Unmapped_filt.R1.fastq
        if [[ -s Unmapped_unfilt.R1.fastq ]] 
        then 
          rm Unmapped_unfilt.*.fastq
        fi
    fi
fi
