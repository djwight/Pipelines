ARG BWAv="0.7.18"
ARG FQCv="0.12.1"
ARG FPv="0.24.0"
ARG SAMv="1.21"
ARG GATKv="4.6.1.0"

FROM buildpack-deps:bookworm AS builder
ARG BWAv FQCv FPv SAMv GATKv

WORKDIR /build

RUN apt-get update && apt-get install -y \
    wget \
    default-jre \
    gcc \
    make \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev

# build BWA
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    git checkout tags/v${BWAv} -b build && \
    make

# get FastQC binary
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FQCv}.zip && \
    unzip fastqc_v${FQCv}.zip && \
    rm fastqc_v${FQCv}.zip

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

# get GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATKv}/gatk-${GATKv}.zip && \
    unzip gatk-${GATKv}.zip && \
    mv gatk-${GATKv} gatk && \
    rm gatk-${GATKv}.zip

# final image
FROM python:3.12 AS final
ARG BWAv FQCv FPv SAMv GATKv
ENV BWAv=$BWAv FQCv=$FQCv FPv=$FPv SAMv=$SAMv GATKv=$GATKv

COPY --from=builder /build /tools

RUN apt-get update && apt-get install -y \
    default-jre \
    r-base

# add tools to PATH
RUN export PATH=/tools:${PATH}

WORKDIR /src

COPY . . 
RUN pip install -r requirements.txt

ENTRYPOINT [ "python3", "pe_fastq_to_cram.py" ]