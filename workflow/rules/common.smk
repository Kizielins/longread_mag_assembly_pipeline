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

