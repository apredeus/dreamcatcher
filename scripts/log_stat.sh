#!/bin/bash 

FLAG=$1

if [[ $FLAG == "preprocessed" ]]
then
    RTOT=`perl -ne 'print "$1" if (m/(\d+) total reads, (\d+) reads not mapped to host genome/)' preprocess_reads.log`
    RUNM=`perl -ne 'print "$2" if (m/(\d+) total reads, (\d+) reads not mapped to host genome/)' preprocess_reads.log`
    RFLT=`perl -ne 'print "$1" if (m/outputted (\d+) reads/)' preprocess_reads.log`
    RDIS=`perl -ne 'print "$1" if (m/discarded (\d+) reads/)' preprocess_reads.log`
    echo -e "\tlog_stat.sh preprocessed: a total of $RTOT reads was mapped to host genome, and $RUNM remained unmapped"
    echo -e "\tlog_stat.sh preprocessed: of these, $RFLT reads will be used, and $RDIS discarded because of homopolymers or length"
    WL=`perl -ne 'print "$1" if (m/found whitelist (.*?);/)' preprocess_reads.log`
    if [[ $WL != "" ]] 
    then
        echo -e "\tlog_stat.sh preprocessed: barcode correction was performed using whitelist $WL"
        RFL2=`perl -ne 'print "$1" if (m/processed fastq file with (\d+) total reads; (\d+) matching a whitelist, (\d+) corrected, and (\d+) unable to correct/)' preprocess_reads.log`
        RWL=`perl -ne 'print "$2" if (m/processed fastq file with (\d+) total reads; (\d+) matching a whitelist, (\d+) corrected, and (\d+) unable to correct/)' preprocess_reads.log`
        RCOR=`perl -ne 'print "$3" if (m/processed fastq file with (\d+) total reads; (\d+) matching a whitelist, (\d+) corrected, and (\d+) unable to correct/)' preprocess_reads.log`
        RUNC=`perl -ne 'print "$4" if (m/processed fastq file with (\d+) total reads; (\d+) matching a whitelist, (\d+) corrected, and (\d+) unable to correct/)' preprocess_reads.log`
        echo -e "\tlog_stat.sh preprocessed: of $RFL2 filtered reads, $RWL matched the whitelist, $RCOR were corrected, and $RUNC could not be corrected"
    else
        echo -e "\tlog_stat.sh preprocessed: no barcodes were found, or no barcode correction needed to be performed"
    fi
elif [[ $FLAG == "detected" ]]
then
    RFLT=`grep "reads; of these"            detected.hisat2_map.log | awk '{print $1}'`
    RKUQ=`grep "superkingdom.*Bacteria$"    kuniq.report.txt        | cut -f2`
    RD1=`grep "aligned exactly 1 time"      detected.hisat2_map.log | awk '{print $1}'`
    RD2=`grep "aligned >1 times"            detected.hisat2_map.log | awk '{print $1}'`
    RD3=$((RD1+RD2))
    ACK=`awk '/superkingdom\t    Bacteria/,/superkingdom\t    Archaea/'  kuniq.report.txt | grep -w assembly | wc -l`
    ACC=`cat detected_strains.list | wc -l`
    ACA=`cut -f12 detected.annotated_fcounts.tsv | sort | uniq | wc -l` ## at least 1 gene detected
    FCD=`awk '{sum+=$7} END {printf "%d",sum}' detected.annotated_fcounts.tsv` ## approximate number of "useful" reads
    GD3=`cut -f1 detected.annotated_fcounts.tsv | wc -l` ## number of detected bacterial genes
    echo -e "\tlog_stat.sh detected: $ACK accessions in KrakenUniq; $ACC detected strains, $ACA detected with at least 1 annotated gene"
    echo -e "\tlog_stat.sh detected: out of $RFLT filtered reads unmapped to host, $RKUQ assigned to Bacteria superkingdom by KrakenUniq"
    echo -e "\tlog_stat.sh detected: $RD3 reads mapped by Hisat2 to all detected strains ($RD1 unique, $RD2 multi-mappers)"
    echo -e "\tlog_stat.sh detected: of these, approximately $FCD reads were assigned to $GD3 bacterial genes"
elif [[ $FLAG == "filtered" ]]
then
    ACF=`cat filtered_strains.list | wc -l`
    FCF=`awk '{sum+=$7} END {printf "%d",sum}' filtered.annotated_fcounts.tsv` ## approximate number of "useful" reads
    GF3=`cut -f1 filtered.annotated_fcounts.tsv | wc -l` ## number of filtered bacterial genes
    NF1=`awk '$11 > 0.1'  detected.annotated_fcounts.tsv | wc -l`
    NF2=`awk '$11 <= 0.1' detected.annotated_fcounts.tsv | awk '($9 == "rRNA" && $10 > 3.4) || ($9 != "rRNA" && $10 > 2.3)' | wc -l` 
    NF3=`awk '$11 <= 0.1' detected.annotated_fcounts.tsv | awk '($9 == "rRNA" && $10 <= 3.4) || ($9 != "rRNA" && $10 <= 2.3)' | grep -cf blacklisted_genes.list`
    echo -e "\tlog_stat.sh filtered: $ACF bacterial strains retained after filtering"
    echo -e "\tlog_stat.sh filtered: $NF1 detected bacterial genes were putatively human, $NF2 had too many mismatches, and $NF3 belonged to the black list" 
    echo -e "\tlog_stat.sh filtered: approximately $FCF reads were assigned to $GF3 filtered bacterial genes."
elif [[ $FLAG == "top" ]]
then
    RFLT=`grep "reads; of these"            top.hisat2_map.log | awk '{print $1}'`
    RT1=`grep "aligned exactly 1 time"      top.hisat2_map.log | awk '{print $1}'`
    RT2=`grep "aligned >1 times"            top.hisat2_map.log | awk '{print $1}'`
    RT3=$((RT1+RT2))
    ACT=`cat top_strains.list | wc -l`
    FCT=`awk '{sum+=$7} END {printf "%d",sum}' top.annotated_fcounts.tsv` ## approximate number of "useful" reads
    GT3=`cut -f1 top.annotated_fcounts.tsv | wc -l` ## number of detected bacterial genes
    echo -e "\tlog_stat.sh top: $ACT bacterial strains retained after the graph collapsing"
    echo -e "\tlog_stat.sh top: out of $RFLT filtered reads unmapped to host, $RT3 reads mapped by Hisat2 to all top strains ($RT1 unique, $RT2 multi-mappers)"
    echo -e "\tlog_stat.sh top: approximately $FCT reads were assigned to $GT3 bacterial genes"
elif [[ $FLAG == "umi" ]] 
then
    GTU=`zcat gene.umi_counts.tsv.gz   | grep -v count | cut -f1 | tr ',' '\n' | sort | uniq | wc -l`
    ACU=`zcat strain.umi_counts.tsv.gz | grep -v count | cut -f1 | sort | uniq | wc -l`
    BCS=`zcat strain.umi_counts.tsv.gz | grep -v count | cut -f2 | sort | uniq | wc -l`
    UMI=`zcat strain.umi_counts.tsv.gz | grep -v count | awk '{sum+=$3} END {printf "%d",sum}'`
    echo -e "\tlog_stat.sh top: $ACU bacterial strains with at least 1 UMI in the final output"
    echo -e "\tlog_stat.sh top: $GTU unique genes spanning $BCS unique barcodes were detected"
    echo -e "\tlog_stat.sh top: a total of $UMI bacterial UMIs detected in all droplets"
else 
    >&2 echo "log_stat.sh: ERROR: $FLAG is not a recognised flag!" 
    exit 1
fi 
