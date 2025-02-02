import os
import logging
from time import time
import tomllib
import json

from modules.logger import init_logger
from modules.utils import run_cmd, validate_file, return_nice_time

# read the config
with open(os.path.join("/fastqs/", os.environ.get("RUN_CONFIG")), "rb") as file:
    cfg = tomllib.load(file)
SAMPLE = cfg["sample"]["name"]

# tool versions from the docker build
VERSIONS = {
    "Minimap2": os.environ.get("MM2v"),
    "FasqQC": os.environ.get("FQCv"),
    "Fastplong": os.environ.get("FPLNGv"),
    "Samtools": os.environ.get("SAMv"),
    "GATK": os.environ.get("GATKv"),
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
    validate_file((fq_path := os.path.join("/fastqs", cfg["sample"]["FQ"])))

    # pre-trimming fastqc report generation
    logging.info(f"Running FastQC...")
    begin = time()
    fastqc_cmd = (
        f"/tools/fastqc --threads={cfg['config']['threads']} -o {OUTDIR} {fq_path}"
    )
    run_cmd(
        cmd=fastqc_cmd,
        tool="fastqc",
        log_stdout=True,
    )
    logging.info(f"FastQC completed in {return_nice_time(begin)}s")

    # read trimming and alignment to genome
    trim_align_cmd = (
        f"/tools/fastplong -w {cfg['config']['threads']} --stdout -h {OUTDIR}/{SAMPLE}.html -j {OUTDIR}/{SAMPLE}.json "
        f"-l {cfg['fastplong']['min_len']} -m {cfg['fastplong']['mean_qual']} -i {fq_path}  "
        f"| /tools/minimap2 -ax map-ont -t {cfg['config']['threads']} /genome/{cfg['ref']['index']} - "
        f"| /tools/st/samtools sort -@ {cfg['config']['threads']} -O bam -T /tmp -o {OUTDIR}/{SAMPLE}.sorted.bam -"
    )

    logging.info(f"Running alignment...")
    begin = time()
    run_cmd(
        cmd=trim_align_cmd,
        tool="fastplong|minimap2|samtools sort",
    )
    logging.info(f"Alignment completed in {return_nice_time(begin, mins=True)}mins")

    # TODO: parse and organise the stats and logging


if __name__ == "__main__":
    runner()
