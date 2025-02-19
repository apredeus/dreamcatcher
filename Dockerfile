# Dockerfile for Seurat 4.3.0
FROM rocker/r-ver:4.2.0

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE
ENV DEBIAN_FRONTEND=noninteractive

# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    libhdf5-dev libcurl4-openssl-dev libssl-dev libpng-dev libboost-all-dev libxml2-dev \
    openjdk-8-jdk python3-dev python3-pip wget git libfftw3-dev libgsl-dev pkg-config \
    pigz zlib1g-dev libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev \
    libssl-dev libdeflate-dev libfontconfig1-dev pbzip2 pigz llvm-10 libgeos-dev 

# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy
RUN pip3 install umap-learn

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

# Install bioconductor dependencies & suggests
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma', 'glmGamPoi'))"

# Install CRAN suggests
RUN R --no-echo --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools'))"

# Install spatstat
RUN R --no-echo --no-restore --no-save -e "install.packages(c('spatstat.explore', 'spatstat.geom'))"

# Install hdf5r
RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"

# Install latest Matrix
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz', repos=NULL, type='source')"

# Install rgeos
RUN R --no-echo --no-restore --no-save -e "install.packages('rgeos')"

# Install Seurat
RUN R --no-restore --no-save -e "install.packages('Seurat')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"

## Install dplyr, igraph 
RUN R --no-echo --no-restore --no-save -e "install.packages(c('dplyr', 'igraph'))"

CMD [ "R" ]
## now install all of the binary tools we need
 
ARG hisat_version=2.2.1
ARG bowtie_version=2.5.1
ARG samtools_version=1.21
ARG umitools_version=1.1.6
ARG kuniq_version=1.0.4
ARG subread_version=2.0.2
ARG bbmap_version=39.15 

#Install hisat2
RUN wget --no-check-certificate https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v${hisat_version}.tar.gz && \
    tar -xvf v${hisat_version}.tar.gz -C /opt && \
    cd /opt/hisat2-${hisat_version} && \
    make && \
    cd / && rm v${hisat_version}.tar.gz

#Install bowtie2
RUN wget --no-check-certificate https://github.com/BenLangmead/bowtie2/archive/refs/tags/v${bowtie_version}.tar.gz && \
    tar -xvf v${bowtie_version}.tar.gz -C /opt && \
    cd /opt/bowtie2-${bowtie_version} && \
    make && \
    cd / && rm v${bowtie_version}.tar.gz

#Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 && \
    tar -xvf samtools-${samtools_version}.tar.bz2 -C /opt && \
    cd /opt/samtools-${samtools_version} && \
    ./configure && \
    make && \
    make install && \
    cd / && rm samtools-${samtools_version}.tar.bz2   

#Install UMItools 
RUN pip3 install umi_tools==$umitools_version && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/local/bin/python

#Install KrakenUniq 
RUN wget --no-check-certificate https://github.com/fbreitwieser/krakenuniq/archive/refs/tags/v${kuniq_version}.tar.gz && \
    tar -xvf v${kuniq_version}.tar.gz -C /opt && \
    cd /opt/krakenuniq-${kuniq_version} && \
    ./install_krakenuniq.sh /usr/local/bin && \
    cd / && rm v${kuniq_version}.tar.gz

#Install edited subread/featureCounts 
RUN cd /opt && \
    git clone https://github.com/apredeus/subread_precision && \
    cd subread_precision/src && \
    make -f Makefile.Linux

## BBMap (for reformat.sh/bbduk.sh)
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_${bbmap_version}.tar.gz && \
    tar -xzf BBMap_${bbmap_version}.tar.gz -C /opt && \
    cd /opt/bbmap && \
    ./stats.sh in=resources/phix174_ill.ref.fa.gz && \
    cd / && rm BBMap_${bbmap_version}.tar.gz

ENV PATH="${PATH}:/opt/hisat2-${hisat_version}:/opt/bowtie2-${bowtie_version}:/opt/subread_precision/bin:/opt/bbmap"     

#Saving Software Versions to a file
RUN echo "hisat2 version: ${hisat_version}" >> versions.txt && \
    echo "bowtie2 version: ${bowtie_version}" >> versions.txt && \
    echo "samtools version: ${samtools_version}" >> versions.txt && \
    echo "UMI-tools version: ${umitools_version}" >> versions.txt && \
    echo "krakenuniq version: ${kuniq_version}" >> versions.txt && \
    echo "subread version: ${subread_version}" >> versions.txt && \
    echo "BBMap version: ${bbmap_version}" >> versions.txt

COPY Dockerfile /docker/
RUN chmod -R 755 /docker
