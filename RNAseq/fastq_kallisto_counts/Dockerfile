ARG FQCv="0.12.1"
ARG FPv="0.24.0"
ARG SAMv="1.21"
ARG KALv="0.51.1"

FROM buildpack-deps:bookworm AS builder
ARG FQCv FPv SAMv KALv

WORKDIR /build

RUN apt-get update && apt-get install -y \
    wget \
    default-jre \
    gcc \
    make \
    cmake \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev

# get FastQC binary
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FQCv}.zip && \
    unzip fastqc_v${FQCv}.zip && \
    rm fastqc_v${FQCv}.zip && \
    mv FastQC/fastqc . && \
    rm -r FastQC


# get fastp binary
RUN wget http://opengene.org/fastp/fastp.${FPv} && \
    mv fastp.${FPv} fastp && \
    chmod a+x ./fastp

# build Samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMv}/samtools-${SAMv}.tar.bz2 && \
    tar -vxjf samtools-${SAMv}.tar.bz2 && \
    rm samtools-${SAMv}.tar.bz2 && \
    mv samtools-${SAMv} st && \
    cd st && \
    ./configure --prefix=/build/st && \
    make

# build kallisto
RUN wget https://github.com/pachterlab/kallisto/archive/refs/tags/v${KALv}.zip && \
    unzip v${KALv}.zip && \ 
    rm v${KALv}.zip && \
    cd kallisto-${KALv} && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    mv src/kallisto ../../ && \
    cd ../.. && \
    rm -r kallisto-${KALv}

# final image
FROM python:3.12 AS final
ARG FQCv FPv SAMv KALv
ENV FQCv=$FQCv FPv=$FPv SAMv=$SAMv KALv=$KALv

COPY --from=builder /build /tools

RUN apt-get update && apt-get install -y \
    default-jre

# add tools to PATH
RUN export PATH=/tools:${PATH}

WORKDIR /pipelines

COPY . . 
RUN pip install -r RNAseq/fastq_kallisto_counts/requirements.txt
RUN python3 setup.py install

ENTRYPOINT [ "python3", "RNAseq/fastq_kallisto_counts/pe_gene_counts.py" ]