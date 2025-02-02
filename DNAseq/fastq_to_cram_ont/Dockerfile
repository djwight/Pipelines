ARG MM2v="2.28"
ARG FQCv="0.12.1"
ARG FPLNGv="0.2.2"
ARG SAMv="1.21"
ARG GATKv="4.6.1.0"

FROM buildpack-deps:bookworm AS builder
ARG MM2v FQCv FPLNGv SAMv GATKv

WORKDIR /build

RUN apt-get update && apt-get install -y \
    wget \
    # default-jre \
    gcc \
    make
# libbz2-dev \
# zlib1g-dev \
# libncurses5-dev \
# libncursesw5-dev \
# liblzma-dev

RUN wget https://github.com/lh3/minimap2/releases/download/v${MM2v}/minimap2-${MM2v}_x64-linux.tar.bz2 && \
    tar -jxvf minimap2-${MM2v}_x64-linux.tar.bz2  && \
    rm minimap2-${MM2v}_x64-linux.tar.bz2 && \
    mv minimap2-${MM2v}_x64-linux/minimap2 . && \
    rm -r minimap2-${MM2v}_x64-linux


# get FastQC binary
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FQCv}.zip && \
    unzip fastqc_v${FQCv}.zip && \
    rm fastqc_v${FQCv}.zip && \
    mv FastQC/fastqc . && \
    rm -r FastQC

# get fastp binary
RUN wget http://opengene.org/fastplong/fastplong.${FPLNGv} && \
    mv fastplong.${FPLNGv} fastplong && \
    chmod a+x ./fastplong

# build Samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMv}/samtools-${SAMv}.tar.bz2 && \
    tar -vxjf samtools-${SAMv}.tar.bz2 && \
    rm samtools-${SAMv}.tar.bz2 && \
    mv samtools-${SAMv} st && \
    cd st && \
    ./configure --prefix=/build/st && \
    make

# # get GATK
# RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATKv}/gatk-${GATKv}.zip && \
#     unzip gatk-${GATKv}.zip && \
#     mv gatk-${GATKv} gatk && \
#     rm gatk-${GATKv}.zip

# final image
FROM python:3.12 AS final
ARG MM2v FQCv FPLNGv SAMv GATKv
ENV MM2v=$MM2v FQCv=$FQCv FPLNGv=$FPLNGv SAMv=$SAMv GATKv=$GATKv

COPY --from=builder /build /tools

# RUN apt-get update && apt-get install -y \
#     default-jre \
#     r-base

# add tools to PATH
RUN export PATH=/tools:${PATH}

WORKDIR /pipelines

COPY . . 
RUN pip install -r DNAseq/fastq_to_cram_ont/requirements.txt
RUN python3 setup.py install

ENTRYPOINT [ "top" ]
#ENTRYPOINT [ "python3", "DNAseq/fastq_to_cram/pe_fastq_to_cram.py" ]