# Fastq to CRAM

Pipeline for alignment of paired-end Illumina short read NGS data to a reference genome from a DNA sequencing experiment. Uses the workhorse BWA-MEM for alignment to a linear reference genome.

This pipeline is set up through docker and runs one sample per docker container. To run multiple samples at once, multiple docker containers can be run simultaneously.

## Configuring the Pipeline Run

The pipeline takes parameters from the `run_config.toml` passed to the container at runtime. 

**Note**- this config file should exist in the directory with the FASTQ files, otherwise it will not be available in the running container.

```toml
[config]
threads = <int:max number of threads the pipeline can use>
debug = <bool:if intermediate files should be saved for debugging>

[sample]
FQ1 = <str:name of fastq read 1>
FQ2 = <str:name of fastq read 2>
known_sites = <str:name of vcf/bcf/bed with known sites of polymorphic sites for BQSR>

[ref]
reference = <str:name of reference file, must have idx in same directory for BWA>

[align]
bwa_flags = <str:optional flags to be given to BWA MEM (see docs)>
bwa_read_group = <str:@RG tag to be given to BWA MEM (see docs)>

[recal_bq]
recal = <bool:if BQSR should be performed (requires sample/known_sites)>
```

## Running the Pipeline
1. Build the docker image.
```bash
docker build -t dnapipe .
```
2. Create a `run_config.toml` (see [Configuring the Pipeline Run](#configuring-the-pipeline-run)). Ensure the toml file is in the same directory as the FASTQ files. Name can be different as long as it is specified in the `-e` environment variable passed to the container at runtime.
3. Run the docker container
```bash
docker run \
    -e RUN_CONFIG=/path/to/config.toml(must be located with FASTQs) \
    --mount type=bind,source=/path/to/genome,target=/genome,readonly \
    --mount type=bind,source=/path/to/fastqs,target=/fastqs \
    -it dnapipe
```
