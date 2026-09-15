"""Shared helpers used across the rule files."""

import csv
import os


def load_samples(sheet_path):
    """Parse the TSV pointed at by config["samples_sheet"] into two dicts:
    fastq (sample -> already-basecalled fastq path) and pod5 (sample -> POD5
    input directory, for samples that need Dorado basecalling first -- see
    basecalling.smk)."""
    fastq, pod5 = {}, {}
    with open(sheet_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            has_fastq = bool(row.get("fastq", "").strip())
            has_pod5 = bool(row.get("pod5_dir", "").strip())
            if has_fastq == has_pod5:
                raise ValueError(
                    f"Sample {sample!r} in {sheet_path!r} must set exactly one of "
                    "`fastq` or `pod5_dir`, not both/neither."
                )
            if has_fastq:
                fastq[sample] = row["fastq"]
            else:
                pod5[sample] = row["pod5_dir"]
    if not fastq and not pod5:
        raise ValueError(f"No samples found in {sheet_path!r} - expected a header row plus one row per sample.")
    return fastq, pod5

def expand_path(path):
    """Expand ~ and $VARS in a config-supplied filesystem path."""
    return os.path.expanduser(os.path.expandvars(path))


def cfg_bool(key, default=False):
    """Read a boolean config value safely."""
    return str(config.get(key, default)).strip().lower() in ("true", "1", "yes")

def active_binners():
    """metabat2 + SemiBin2 always; COMEBin only if use_comebin is set """
    binners = ["metabat2", "semibin2"]
    if cfg_bool("use_comebin"):
        binners.append("comebin")
    return binners

def get_bin_ids(wildcards):
    """Dynamic list of DAS_Tool-refined bin IDs (bin.1, bin.2, ...).
    Any rule that needs all bins must call this function.
    """
    refined_dir = checkpoints.das_tool_refine.get(**wildcards).output.bins_dir
    bin_ids = glob_wildcards(os.path.join(refined_dir, "{bin}.fa")).bin
    return sorted(bin_ids, key=lambda b: int(b.split(".")[-1]))
