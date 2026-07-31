"""Shared helpers used across the rule files."""

import csv
import os


def load_samples(sheet_path):
    """Parse the TSV pointed at by config["samples_sheet"] into {sample: fastq}."""
    samples = {}
    with open(sheet_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            samples[row["sample"]] = row["fastq"]
    if not samples:
        raise ValueError(f"No samples found in {sheet_path!r} -- expected a header row plus one row per sample.")
    return samples


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