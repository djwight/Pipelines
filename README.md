# Pipelines

Collection of various bioinformatics pipelines for DNA and RNA sequencing.

## Repository Structure

Pipelines for DNA sequencing and RNA sequencing tasks can be found in the sub-directories highlighted below.

```
.
├── DNAseq
│   ├── fastq_to_cram
│   └── fastq_to_cram_ont
├── modules
├── RNAseq
│   └── fastq_kallisto_counts
└── setup.py
```

## Tests

Unit tests are run using `pytest` and can be run with the following command in the root directory. 

```bash
pytest tests
```

Configuration for the tests can be found in the `.pytest.ini` configuration file. See the [pytest documentation](https://docs.pytest.org/en/stable/reference/customize.html#configuration-file-formats) for more info.