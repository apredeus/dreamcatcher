# Dreamcatcher

Dreamcatcher is an advanced pipeline for accurate identification of _bacterial_ reads in human RNA-seq data. It has two modes: (1) bulk, in which it identifies bacterial species present in the sample and quantifies their abundance & gene expression; (2) single cell, which includes several extra steps which, in addition to the "bulk" output, result in a (barcode) x (bacterial species) matrix, which could be conveniently added to the single cell/spatial metadata. 

Unfortunately, only 10x protocols are currently supported on a single-cell level. You can run Smart-seq(1/2/3) or other plate-based data in a bulk mode though. 

## Motivation

Main motivation behind `dreamcatcher` is as follows: most metagenome profiling methods can be roughtly divided into k-mer based and mapping-based. Mapping-based methods are usually considered to be too slow to run against a modern bacterial database, containing 10's of thousands of genomes. At the same time, k-mer based methods produce overwhelming number of false positives, which come from two sources: (1) spurious matches to foreign sequences in the reference genomes; (2) huge sampling bias towards certain bacterial species that are human pathogens. Both problems are exacerbated in case of RNA-seq, since most (~90%) of the detected reads come from bacterial rRNA operons. 

At the same time, careful case-by-case consideration of bacteria of interest has shown that it's possible to distinguish spurious and proper alignments when using mapped reads and accurate genome annotation. Furthermore, k-mer-based algorithms do not suffer from substantial false negative rates; simply put, if the bug is present in the reference and in the sample, they would detect it. We decided to use a combination of these methods, and further leverage gene annotation to reliably define the bacterial species truly present in the sample. 

## Algorithm

`dreamcatcher` leverages RefSeq genomes with matched curated genome annotations (see section on the Reference preparation below). Our current reference follows regular Kraken/Kraken2/KrakenUniq approach. This allows us to put a comprehensive RefSeq-based 

  - reads are mapped using STAR, STARsolo, Cell Ranger, or a similar tool. 
  - in case of 10x single cell/Visium protocol, unmapped reads are converted into a "UMI-tools" style of formatting, with barcode and UMI sequences added to the read name; 
  - KrakenUniq is used to classify the unmapped reads, and assemblies with at least minimal number of unique k-mers are selected (usually 200-1000 strains); 
  - Hisat2 is used to map the same unmapped reads to this selection of genomes, and featureCounts is used to find the expressed genes; 
  - expression table is then annotated according to the gene type (rRNA or not), mean number of mismatches per 100 bp, etc;
  - following this, unmapped reads are mapped to T2T human assembly as well as PhiX; 
  - genes that have too many mismatches, that map to human assembly or PhiX, or genes from the blacklist are removed; 
  - strains, that do not express a minimal number of genes are removed as contaminants/artifacts; 
  - following this, a graph of shared reads is generated for the filtered strains; 
  - the graph is clustered iteratively, based on shared read topology and taxonomy;
  - this results in selection of strains named _top strains_; this is the main result of the "bulk" part of the pipeline;
  - in case of single-cell mode, the unmapped reads are then mapped to the genomes of _top strains_, and the resulting BAM file is counted using `featureCounts`; 
  - using `UMI-tools` allows us to get UMI counts per barcode, per _top strains_ gene; 
  - finally, we sum the expression up per species and generate a final (barcode) x (species) table. 

## Output 

The output contains the following folders: 

## Reference preparation 

Reference is prepared 
