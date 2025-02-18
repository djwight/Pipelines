from pathlib import Path
from modules.stats import create_dna_run_stats, create_rna_pseudo_stats


def test_dna_stats():
    stats = create_dna_run_stats(
        tool_versions={"tool1": "3.4.1", "tool2": "0.1.3"},
        sample="sample",
        outdir=Path(__file__).parent / "test_data",
    )
    assert all(
        [i in ["tools", "fastq", "alignment", "deduplication"] for i in stats.keys()]
    )
    assert stats["tools"]["tool1"] == "3.4.1"
    assert "fastp_version" in stats["fastq"].keys()
    assert "filtered_sequences" in stats["alignment"].keys()
    assert "READ_PAIRS_EXAMINED" in stats["deduplication"].keys()


def test_rna_stats():
    stats = create_rna_pseudo_stats(
        tool_versions={"tool1": "3.4.1", "tool2": "0.1.3"},
        sample="sample",
        outdir=Path(__file__).parent / "test_data",
    )
    assert all([i in ["tools", "fastq", "pseudo_alignment"] for i in stats.keys()])
    assert stats["tools"]["tool1"] == "3.4.1"
    assert "fastp_version" in stats["fastq"].keys()
    assert "kallisto_version" in stats["pseudo_alignment"].keys()
