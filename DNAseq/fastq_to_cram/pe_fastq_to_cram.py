import os
import logging
from time import time
import re
import tomllib
import json

from modules.logger import init_logger
from modules.utils import run_cmd, validate_file, return_nice_time
from modules.stats import create_dna_run_stats

# read the config
with open(os.path.join("/fastqs/", os.environ.get("RUN_CONFIG")), "rb") as file:
    cfg = tomllib.load(file)
SAMPLE = re.split(r"_\d.", cfg["sample"]["FQ1"])[0]

# tool versions from the docker build
VERSIONS = {
    "BWA": os.environ.get("BWAv"),
    "FasqQC": os.environ.get("FQCv"),
    "Fastp": os.environ.get("FPv"),
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
    validate_file((fq1_path := os.path.join("/fastqs", cfg["sample"]["FQ1"])))
    validate_file((fq2_path := os.path.join("/fastqs", cfg["sample"]["FQ2"])))

    # pre-trimming fastqc report generation
    logging.info(f"Running FastQC...")
    begin = time()
    fastqc_cmd = f"/tools/FastQC/fastqc --threads={cfg['config']['threads']} -o {OUTDIR} {fq1_path} {fq2_path}"
    run_cmd(
        cmd=fastqc_cmd,
        tool="fastqc",
        log_stdout=True,
    )
    logging.info(f"FastQC completed in {return_nice_time(begin)}s")

    # read trimming and alignment to genome
    trim_align_cmd = (
        f"/tools/fastp -w {cfg['config']['threads']} --stdout -h {OUTDIR}/{SAMPLE}.html -j {OUTDIR}/{SAMPLE}.json -i {fq1_path} -I {fq2_path} "
        f"| /tools/bwa/bwa mem {cfg['align']['bwa_options']} -M -p -t {cfg['config']['threads']} -R '{cfg['align']['bwa_read_group']}' "
        f"/genome/{cfg['ref']['reference']} - | /tools/st/samtools sort -@ {cfg['config']['threads']} -O bam -T /tmp -o {OUTDIR}/{SAMPLE}.sorted.bam -"
    )

    logging.info(f"Running alignment...")
    begin = time()
    run_cmd(
        cmd=trim_align_cmd,
        tool="fastp|bwa mem|samtools sort",
    )
    logging.info(f"Alignment completed in {return_nice_time(begin, mins=True)}mins")

    # Gather alignment stats
    align_stats_cmd = f"/tools/st/samtools stats {OUTDIR}/{SAMPLE}.sorted.bam | grep ^SN | cut -f 2- > {OUTDIR}/{SAMPLE}.aln.stats"

    logging.info(f"Gathering alignment stats...")
    begin = time()
    run_cmd(
        cmd=align_stats_cmd,
        tool="samtools stats",
    )
    logging.info(
        f"Stats generation completed in {return_nice_time(begin, mins=True)}mins"
    )

    # mark duplicates
    logging.info(f"Removing duplicates...")
    rm_dups_cmd = (
        f"/tools/gatk/gatk MarkDuplicates -R /genome/{cfg['ref']['reference']} -I {OUTDIR}/{SAMPLE}.sorted.bam "
        f"-O {OUTDIR}/{SAMPLE}.sorted.dedup.bam -M {OUTDIR}/{SAMPLE}_dup_metrics.txt"
    )
    begin = time()
    run_cmd(cmd=rm_dups_cmd, tool="MarkDuplicates")
    logging.info(f"Duplicates removed in {return_nice_time(begin, mins=True)}s")

    # recalibrate BQs
    if cfg["recal_bq"]["recal"]:
        logging.info(f"Recalibrating BQ scores...")
        bqsr_cmd = (
            f"/tools/gatk/gatk BaseRecalibrator -R /genome/{cfg['ref']['reference']} --known-sites {cfg['sample']['known_sites']} "
            f"-I {OUTDIR}/{SAMPLE}.sorted.dedup.bam -O {OUTDIR}/{SAMPLE}_recal_data.table"
        )
        apply_bqsr_cmd = (
            f"/tools/gatk/gatk ApplyBQSR -R /genome/{cfg['ref']['reference']} -I {OUTDIR}/{SAMPLE}.sorted.dedup.bam "
            f"--bqsr-recal-file {OUTDIR}/{SAMPLE}_recal_data.table -O {OUTDIR}/{SAMPLE}.recal.sorted.dedup.bam && "
            f"mv {OUTDIR}/{SAMPLE}.recal.sorted.dedup.bam {OUTDIR}/{SAMPLE}.sorted.dedup.bam"
        )
        begin = time()
        run_cmd(cmd=bqsr_cmd, tool="BaseRecalibrator", log_stdout=True)
        run_cmd(cmd=apply_bqsr_cmd, tool="ApplyBQSR", log_stdout=True)
        logging.info(f"BQSR applied in {return_nice_time(begin, mins=True)}mins")

    # convert BAM to CRAM
    logging.info(f"Converting BAM to CRAM...")
    cram_cmd = (
        f"/tools/st/samtools view -@ {cfg['config']['threads']} --cram -T /genome/{cfg['ref']['reference']} "
        f"-o {OUTDIR}/{SAMPLE}.sorted.dedup.cram {OUTDIR}/{SAMPLE}.sorted.dedup.bam"
    )
    begin = time()
    run_cmd(cmd=cram_cmd, tool="samtools view", log_stdout=True)
    logging.info(f"CRAM created in {return_nice_time(begin, mins=True)}mins")

    # write stats file
    stats = create_dna_run_stats(tool_versions=VERSIONS, sample=SAMPLE, outdir=OUTDIR)
    with open(f"{OUTDIR}/{SAMPLE}.stats", "w") as file:
        json.dump(stats, file)

    # clean up
    if not cfg["config"]["debug"]:
        os.remove(f"{OUTDIR}/{SAMPLE}.sorted.dedup.bam")
        os.remove(f"{OUTDIR}/{SAMPLE}.sorted.bam")

    logging.info(f"Script completed in {return_nice_time(start, mins=True)}mins")


if __name__ == "__main__":
    runner()
