"""Filter bins by CheckM2 completeness/contamination stats and copy 
to Final-bins/passed/ for GTDB-Tk classification. """

import csv
import os
import shutil

min_completeness = float(snakemake.params.min_completeness)
max_contamination = float(snakemake.params.max_contamination)

os.makedirs(snakemake.output.finalbins, exist_ok=True)

bins_by_id = {os.path.splitext(os.path.basename(p))[0]: p for p in snakemake.input.bins}

passed = []
with open(snakemake.input.checkm2_report) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    for row in reader:
        if len(row) < 3:
            continue
        bin_id, completeness, contamination = row[0], float(row[1]), float(row[2])
        if completeness >= min_completeness and contamination <= max_contamination:
            passed.append(bin_id)
            shutil.copy(bins_by_id[bin_id], os.path.join(snakemake.output.finalbins, f"{bin_id}.fasta"))

with open(snakemake.output.passed_list, "w") as fh:
    fh.write("\n".join(passed) + ("\n" if passed else ""))

if not passed:
    raise SystemExit(
        "No bins passed the quality filter "
        f"(completeness >= {min_completeness}%, contamination <= {max_contamination}%)."
    )
