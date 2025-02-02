# RNAseq: Fastq to Transcript Counts (Pair-ended Sequencing)

Pipeline for quantification of transcript counts from pair-end sequencing data. Uses the pseudo-alignment algorithm in Kallisto to establish which transcript the reads belong to in the transcriptome provided for the sequenced organism.

This pipeline was designed with Illumnina PE sequencing data in mind, generated from a transcriptome sequencing-appropriate wet-lab experiment. 

This pipeline is set up through docker and runs one sample per docker container. To run multiple samples at once, multiple docker containers can be run simultaneously.

## Configuring the Pipeline Run

The pipeline takes parameters from a `run_config.toml` passed to the container at runtime. See `example_config.toml` and below to understand the options available. 

**Note**- this config file should exist in the directory with the FASTQ files, otherwise it will not be available in the running docker container.

```toml
[config]
threads = <int:max number of threads the pipeline can use>
debug = <bool:if intermediate files should be saved for debugging>

[sample]
name = <str:name of the sample>
FQ1 = <str:name of fastq read 1>
FQ2 = <str:name of fastq read 2>

[ref]
index = <str:name of index file for the transcriptome file>
```

## Inputs

This pipeline requires the following inputs to run:

1. 2x FASTQ files (mandatory): from the PE RNA sequencing experiment (names in `run_config.toml`).
2. A transcriptome fasta (mandatory)
3. An index for the transcriptome fasta (mandatory): corresponding to (2) made using Kallisto (see below for more info; name in `run_config.toml`).


### Making an Index File for the Transcciptome Fasta

An index file for the transctiptome fasta can be made using the command below (using the default 31 length kmers).
```bash
kallisto index -t <number_of_threads> -i <name_of_index_file> <path_to_transcriptome_fasta>
```

For more info see [Kallisto docs](https://pachterlab.github.io/kallisto/manual).


## Running the Pipeline
1. Build the docker image from upper directory (docker context needs access to the modules in `src/`).
```bash
docker build -t kalpipe -f RNAseq/fastq_kallisto_counts/Dockerfile .
```
2. Create a `run_config.toml` (see [Configuring the Pipeline Run](#configuring-the-pipeline-run)). Ensure the toml file is in the same directory as the FASTQ files. Name can be different as long as it is specified in the `-e` environment variable passed to the container at runtime.
3. Run the docker container
```bash
docker run \
    -e RUN_CONFIG=<name of config.toml (must be located with FASTQs)> \
    --mount type=bind,source=/path/to/genome,target=/genome,readonly \
    --mount type=bind,source=/path/to/fastqs,target=/fastqs \
    -it kalpipe
```

## Outputs

A directory will be created in the directory for the FASTQ files called `<sample_name>_result/` where the results will be posted. Results directory is organized as below and files are described.

```
├─ <sample_name>_result
   ├── abundance.tsv => output from kallisto quant.
   ├── <sample_name>.html => output from fastp as html.
   ├── <sample_name>.json => output from fastp as json.
   ├── <sample_name>.stats => custom stats collected from the experiment.
   ├── run_info.json => stats output by jllisto quant.
   └── run.log => logs capture from the docker container.
```


**Note-** As docker runs containers as root by standard, the following will need to be done to change the owner to the current user (on linux systems):

```bash
sudo chown -R user:user <sample_name>_result/
```