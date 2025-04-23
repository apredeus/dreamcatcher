# Dreamcatcher

<img src="https://github.com/apredeus/dreamcatcher/blob/main/img/dreamcatcher.png">

Dreamcatcher is an advanced pipeline for accurate identification of _bacterial_ reads in human RNA-seq data. It has two modes: 

  - bulk, in which it identifies bacterial species present in the sample and quantifies their abundance & gene expression; 
  - single cell, which includes several extra steps which, in addition to the "bulk" output, result in a (barcode) x (bacterial species) matrix, which could be conveniently added to the single cell/spatial metadata. 

Unfortunately, only 10x protocols are currently supported on a single-cell level. You can run Smart-seq(1/2/3) or other plate-based data in a bulk mode though. 

## Motivation

Main motivation behind `dreamcatcher` is as follows: most metagenome profiling methods can be roughtly divided into k-mer based and mapping-based. Mapping-based methods are usually considered to be too slow to run against a modern bacterial database, containing 10's of thousands of genomes. At the same time, k-mer based methods produce overwhelming number of false positives, which come from two sources: (1) spurious matches to foreign sequences in the reference genomes; (2) huge sampling bias towards certain bacterial species that are human pathogens. Both problems are exacerbated in case of RNA-seq, since most (~90%) of the detected reads come from bacterial rRNA operons. 

At the same time, careful case-by-case consideration of bacteria of interest has shown that it's possible to distinguish spurious and proper alignments when using mapped reads and accurate genome annotation. Furthermore, when only minimal filtering is applied, k-mer-based algorithms do not suffer from substantial false negative rates; simply put, if the bug is present in the reference and in the sample, they would detect it. We decided to use a combination of these methods, and further leverage gene annotation to reliably define the bacterial species truly present in the sample. 

## Algorithm

`dreamcatcher` leverages RefSeq genomes with matched curated genome annotations (see section on the Reference preparation below). Our current reference follows regular Kraken/Kraken2/KrakenUniq approach. This allows us to put a comprehensive RefSeq-based 

  - reads are mapped using STAR, STARsolo, Cell Ranger, or a similar tool. 
  - in case of 10x single cell/Visium protocol, unmapped reads are converted into a "UMI-tools" style of formatting, with barcode and UMI sequences added to the read name; 
  - KrakenUniq is used to classify the unmapped reads, and assemblies with at least minimal number of unique k-mers are selected (usually 200-1000 strains). These are the _detected strains_;
  - Hisat2 is used to map the same unmapped reads to this selection of genomes, and featureCounts is used to find the expressed genes; 
  - expression table is then annotated according to the gene type (rRNA or not), mean number of mismatches per 100 bp, etc;
  - following this, unmapped reads are mapped to T2T human assembly as well as PhiX; 
  - genes that have too many mismatches, that map to human assembly or PhiX, or genes from the blacklist are removed; 
  - strains which do not express a minimal number of genes are removed as contaminants/artifacts. What remains is referred to as _filtered strains_;
  - following this, a graph of shared reads is generated for the filtered strains; 
  - the graph is clustered iteratively, based on shared read topology and taxonomy;
  - this results in selection of strains named _top strains_; this is the main result of the "bulk" part of the pipeline;
  - in case of single-cell mode, the unmapped reads are then mapped to the genomes of _top strains_, and the resulting BAM file is counted using `featureCounts`; 
  - using `UMI-tools` allows us to get UMI counts per barcode, per _top strains_ gene; 
  - finally, we sum the expression up per species and generate a final (barcode) x (species) table. 

## Output 

The output contains the following folders: 

  - `0_preprocessed_reads` - with read preprocessing logs and barcode corrections (if performed); 
  - `1_krakenuniq` - the output of `krakenuniq` unmapped read classification; 
  - `2_detected_strains` - list of the detected strains, as well as `hisat2` indexes, BAM files, annotated `featureCounts` outputs, and per-read annotations;
  - `3_filtered_strains` - list of the filtered strains, as well as annotated per-gene and per-read outputs; 
  - `4_read_networks` - all the files that are the input and the output of the read graph collapsing algorithm; 
  - `5_top_strains` - list of the top strains, as well as `hisat2` indexes, BAM files, annotated `featureCounts` outputs; 
  - `6_umitools` - results of `umi-tools` processing, in case of 10x single cell or Visium experiments. Not present if the tool was ran in the `--bulk` mode.  

Most important output files are also present in the root directory:
  
  - `kuniq.report.txt` - KrakenUniq report file that contains metagenomic classification for all reads unmapped to host, and can be visualised using [Pavian](https://github.com/fbreitwieser/pavian); 
  - `detected.annotated_fcounts.tsv`, `detected.summary.tsv` - per-gene and per-strain metrics for the detected strains, derived by mapping reads to a set of strains (`2_detected_strains/detected_strains.list`) with > 3 unique k-mers in KrakenUniq;  
  - `filtered.annotated_fcounts.tsv`, `filtered.summary.tsv` - per-gene and per-strain metrics for the filtered strains. Details of the filtering process could be found in `3_filtered_strains/filter_strains.log`; 
  - `filtered.cluster.tsv`, `top.cluster.tsv` - results of read network clustering (input and output strains);
  - `mapstats.tsv` - detailed statistics of reads that were mapped/assigned to a gene in detected, filtered, and top strains; 
  - `barcode_to_species.tsv` - corrected barcode by identified bacterial species matrix. Not present if the tool was ran in the `--bulk` mode.

Two file types that don't have a header are `*.summary.tsv` and `*.annotated_fcounts.tsv`. Given below is the description of each file. 

Files `detected.annotated_fcounts.tsv`, `filtered.annotated_fcounts.tsv`, and `top.annotated_fcounts.tsv` all are made using `featureCounts` output, and extra annotation for each gene was added to them. There are 14 columns overall, and following table explains each one using an example output: 

<div align="center">

| RefSeq locus tag | Contig | Gene start | Gene end | Strand | Gene length | Split read count | Raw read count | Gene biotype | Mismatch rate | Fraction T2T-human/phiX | RefSeq strain ID | Strain name | Strain taxid | 
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
P9N54_RS23305 | NZ_CP121294.1 | 4776927 | 4778102 | - | 1176 | 0.142853 | 1 | protein_coding | 0.662 | 0.000 | GCF_029674705.1 | Escherichia coli O155:H21 strain=NWU_1 | 3038394 |
P9N54_RS23360 | NZ_CP121294.1 | 4788469 | 4790003 | - | 1535 | 20.634720 | 6231 | rRNA | 0.516 | 0.000 | GCF_029674705.1 | Escherichia coli O155:H21 strain=NWU_1 | 3038394 |
NX106_RS23395 | NZ_JAMSJK010000099.1 | 1 | 321 | + | 321 | 0.125000 | 1 | rRNA | 2.532 | 1.000 | GCF_024733345.1 | Escherichia ruysiae strain=C61-1 | 2608867 |
NX106_RS23520 | NZ_JAMSJK010000126.1 | 194 | 280 | + | 87 | 0.003418 | 4 | rRNA | 0.625 | 0.000 | GCF_024733345.1 | Escherichia ruysiae strain=C61-1 | 2608867 |

</div>

## Reference preparation 

Dreamcatcher reference takes a while to make. The list of used RefSeq strains, together with the assembly metadata, can be found in this repository (`/data/assembly_summary_bacteria.txt.gz`). Two of the biggest items in the referense are (1) download of RefSeq genomes and annotations; (2) creation of the KrakenUniq database; (3) creation of T2T+phiX decoy `hisat2` index. 

To get RefSeq genomes and annotations in GFF3 format, you would need to install [NCBI datasets and dataformat utilities](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/). Following this, take the strain list (`/data/refseq_strain.list` in this repo) and run

```bash
datasets download genome accession --include genome,gff3,gtf --dehydrated --inputfile refseq_strain.list
unzip ncbi_dataset.zip
datasets rehydrate --gzip  --directory .
```
This will create a directory with subdirectories named `GCF_007970465.1` etc; each individual subdirectory should contain the gzipped genome file (ending with `genomic.fna.gz`) and GFF3 file (`genomic.gff.gz`). Make sure all strains have these two files downloaded successfully before moving on. The current database has 51,800 bacterial species, and takes up about 200 Gb of space.  

To create `krakenuniq` database, we have used the following commands. To re-create the exact database used in the paper, the strain list above needs to be used. The command below used 32 cores and limited the database t o 480 Gb for sensitivity. The final folder needs to be cleaned up and takes up about 600 Gb. 

```bash
krakenuniq-download --threads 32 --db . 'refseq/vertebrate_mammalian/Any/species_taxid=9606'
krakenuniq-download --threads 32 --db . refseq/viral/Any viral-neighbors
krakenuniq-download --threads 32 --db . --dust refseq/archaea
krakenuniq-download --threads 32 --db . --dust refseq/bacteria/Contig/refseq_category='representative genome'
krakenuniq-download --threads 32 --db . --dust refseq/bacteria/Contig/refseq_category='reference genome'
krakenuniq-download --threads 32 --db . --dust refseq/bacteria/Scaffold/refseq_category='representative genome'
krakenuniq-download --threads 32 --db . --dust refseq/bacteria/Scaffold/refseq_category='reference genome'

krakenuniq-build --jellyfish-bin /software/cellgeni/miniconda3/bin/jellyfish --max-db-size 480 --db . --kmer-len 31 --threads 32 --taxids-for-genomes --taxids-for-sequences
```
This takes a long time to create (approximately 2 weeks on 32 cores). Please create a free [Globus](https://www.globus.org/) account and let me know the email you've registered with if you need to download the pre-made reference. 

Additional files required by `dreamcatcher` are: 

  - `seqid2taxid.map.refseq` - full mapping of all sequences available in our KrakenUniq database. These are generated by `krakenuniq` when making the reference; 
  - `assembly_summary_bacteria.txt` - GenBank metadata of all used RefSeq bacterial assemblies, including assembly level and strain info; 
  - `gcf_species_genus.tsv` - mapping of all bacterial RefSeq IDs to NCBI species/genus names and TaxIDs. 

All of these smaller files are available as archives in `/data` subdirectory of this repo. 

## Pipeline organization 

The pipeline is linear, and written in `bash`, with `perl` used extensively for text file parsing. All of the binary executables are called via `singularity`. Thus, you should be able to run it on any Unix system with `singularity` installed. 

### Installation 

To run the pipeline: 

  - create all of the necessary reference indexes (see #Reference preparation above);
  - run `git clone https://github.com/apredeus/dreamcatcher`, and add the folder location to your `$PATH`; 
  - edit the `dreamcatcher` script to include the relevant reference index locations. 

You should be all set after this! 

### Running the analysis

To run the analysis, simply run 

```bash
dreamcatcher <path_to_fastq_directory> <sample_ID> 
```

The default resources required by the script are 16 cores and 128 GB of RAM. You probably can try to run it on as low as 30 GB, depending on the number of detected strains. 
