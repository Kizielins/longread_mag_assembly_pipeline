"""Stage 2 of bin filtering: drop bins GUNC flags as chimeric (pass.GUNC ==
False) from filter_checkm2_bins.py's completeness/contamination survivors,
and copy the final set into staging/passed-bins/ for GTDB-Tk.
"""

import csv
import glob
import os
import shutil
import sys

with open(snakemake.input.checkm2_passed_list) as fh:
    checkm2_passed = [line.strip() for line in fh if line.strip()]

if not checkm2_passed:
    os.makedirs(snakemake.output.finalbins, exist_ok=True)
    open(snakemake.output.passed_list, "w").close()
    print("No CheckM2-passed bins to GUNC-filter.", file=sys.stderr)
    raise SystemExit(0)

matches = glob.glob(os.path.join(snakemake.input.gunc_report, "*maxCSS_level.tsv"))
if not matches:
    raise SystemExit(f"gunc_run produced no *maxCSS_level.tsv in {snakemake.input.gunc_report!r}.")
gunc_pass = {}
with open(matches[0]) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        bin_id = row["genome"].removesuffix(".fasta")
        gunc_pass[bin_id] = row["pass.GUNC"].strip().lower() == "true"

os.makedirs(snakemake.output.finalbins, exist_ok=True)

passed = []
for bin_id in checkm2_passed:
    if not gunc_pass.get(bin_id, False):
        continue
    passed.append(bin_id)
    shutil.copy(
        os.path.join(snakemake.input.bins, f"{bin_id}.fasta"),
        os.path.join(snakemake.output.finalbins, f"{bin_id}.fasta"),
    )

with open(snakemake.output.passed_list, "w") as fh:
    fh.write("\n".join(passed) + ("\n" if passed else ""))

if not passed:
    print ("No bins passed quality filtering (completeness/contamination and GUNC).", file=sys.stderr)
