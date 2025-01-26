import json


def create_dna_run_stats(tool_versions: dict, sample: str, outdir: str) -> dict:
    """Creates a summary stats dictionary from the stats files created by the tools used

    Args:
        tool_versions (dict): versions of tools used in the analysis
        sample (str): name of the sample
        outdir (str): path to the directory for the results

    Returns:
        dict: stats obtained from the different bioinformatics tools
    """
    stats = {}
    stats["tools"] = tool_versions

    # get fastq stats
    with open(f"{outdir}/{sample}.json", "r") as handle:
        stats["fastq"] = json.load(handle)["summary"]

    # get alignment stats
    with open(f"{outdir}/{sample}.aln.stats") as handle:
        tmp_data = [
            line.rstrip().replace(":", "").split("\t")[:2]
            for line in handle.readlines()
        ]
        stats["alignment"] = {k.replace(" ", "_"): v for k, v in tmp_data}

    # get deduplication stats
    keep = False
    tmp_data = []
    with open(f"{outdir}/{sample}_dup_metrics.txt", "r") as handle:
        for line in handle.readlines():
            if line.startswith("## METRICS CLASS"):
                keep = True
                continue
            elif line.startswith("\n"):
                keep = False
            if keep:
                tmp_data.append(line.rstrip().split("\t"))
    stats["deduplication"] = {k: v for k, v in list(zip(tmp_data[0], tmp_data[1]))}
    return stats


def create_rna_pseudo_stats(tool_versions: dict, sample: str, outdir: str) -> dict:
    stats = {}
    stats["tools"] = tool_versions

    # get fastq stats
    with open(f"{outdir}/{sample}.json", "r") as handle:
        stats["fastq"] = json.load(handle)["summary"]

    # get kallisto stats
    with open(f"{outdir}/run_info.json", "r") as handle:
        stats["pseudo_alignment"] = json.load(handle)
    return stats
