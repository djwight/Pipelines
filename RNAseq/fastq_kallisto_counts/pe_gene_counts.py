import os
import logging
from time import time
import re
import tomllib
import json

from modules.logger import init_logger
from modules.utils import run_cmd, validate_file, return_nice_time
from modules.stats import create_rna_pseudo_stats

# read the config
with open(os.path.join("/fastqs/", os.environ.get("RUN_CONFIG")), "rb") as file:
    cfg = tomllib.load(file)
SAMPLE = re.split(r"_\d.", cfg["sample"]["FQ1"])[0]

# tool versions from the docker build
VERSIONS = {
    "FasqQC": os.environ.get("FQCv"),
    "Fastp": os.environ.get("FPv"),
    "Samtools": os.environ.get("SAMv"),
    "Kallisto": os.environ.get("KALv"),
}
os.makedirs((OUTDIR := f"/fastqs/{SAMPLE}_result"), exist_ok=True)

init_logger(outdir=OUTDIR)


def runner() -> None:
    start = time()
    logging.info(f"Script started for sample = {SAMPLE}")

    # log out tool versions used in script
    for tool, version in VERSIONS.items():
        logging.info(f"{tool} version = {version}")

    logging.info("Checking fastq files exist...")
    validate_file((fq1_path := os.path.join("/fastqs", cfg["sample"]["FQ1"])))
    validate_file((fq2_path := os.path.join("/fastqs", cfg["sample"]["FQ2"])))

    # pre-trimming fastqc report generation
    logging.info(f"Running FastQC...")
    begin = time()
    fastqc_cmd = f"/tools/fastqc --threads={cfg['config']['threads']} -o {OUTDIR} {fq1_path} {fq2_path}"
    run_cmd(
        cmd=fastqc_cmd,
        tool="fastqc",
        log_stdout=True,
    )
    logging.info(f"FastQC completed in {return_nice_time(begin)}s")

    # read trimming
    fq1_trim_path = "_trim.".join(fq1_path.split("."))
    fq2_trim_path = "_trim.".join(fq2_path.split("."))
    trim_cmd = (
        f"/tools/fastp -w {cfg['config']['threads']} --stdout -h {OUTDIR}/{SAMPLE}.html -j {OUTDIR}/{SAMPLE}.json "
        f"-i {fq1_path} -I {fq2_path} -o {fq1_trim_path} -O {fq2_trim_path}"
    )
    logging.info(f"Running trimming...")
    begin = time()
    run_cmd(
        cmd=trim_cmd,
        tool="fastp",
    )
    logging.info(f"Trimming completed in {return_nice_time(begin)}s")

    # pseudoalignment and transcript counting
    pseudo_cmd = f"/tools/kallisto quant -t 8 -i {cfg['ref']['index']} -o {OUTDIR} {fq1_trim_path} {fq2_trim_path}"
    logging.info(f"Running trimming...")
    begin = time()
    run_cmd(
        cmd=pseudo_cmd,
        tool="kallisto quant",
    )
    logging.info(
        f"Pseudoalignment completed in {return_nice_time(begin, mins=True)}mins"
    )

    # write stats file
    stats = create_rna_pseudo_stats(
        tool_versions=VERSIONS, sample=SAMPLE, outdir=OUTDIR
    )
    with open(f"{OUTDIR}/{SAMPLE}.stats", "w") as file:
        json.dump(stats, file)

    logging.info(f"Script completed in {return_nice_time(start, mins=True)}mins")


if __name__ == "__main__":
    runner()
