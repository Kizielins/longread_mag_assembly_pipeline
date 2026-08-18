"""Build pipeline.ont.genome.summary and the published Final-bins/ FASTAs."""


import csv
import glob
import os

os.makedirs(snakemake.output.finalbins, exist_ok=True)

with open(snakemake.input.passed_list) as fh:
    passed_bins = [line.strip() for line in fh if line.strip()]
passed_bins.sort(key=lambda b: int(b.split(".")[-1]))

## CheckM2: bin_id -> {completeness, contamination, genome_size, n_contig, n50, gc}
checkm2 = {}
with open(snakemake.input.checkm2_report) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        checkm2[row["Name"]] = {
            "completeness": row.get("Completeness", "NA"),
            "contamination": row.get("Contamination", "NA"),
            "genome_size": row.get("Genome_Size", "NA"),
            "n_contig": row.get("Total_Contigs", "NA"),
            "n50": row.get("Contig_N50", "NA"),
            "gc": row.get("GC_Content", "NA"),
        }

## GTDB-Tk: bin_id -> taxonomy string
taxa = {}
with open(snakemake.input.taxa) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            taxa[parts[0]] = parts[1]

## contigs_info.tsv (assembly.smk): contig_id -> circular
circular_by_contig = {}
with open(snakemake.input.contigs_info) as fh:
    reader = csv.reader(fh, delimiter="\t")
    next(reader)  # header
    for row in reader:
        if len(row) >= 3:
            circular_by_contig[row[0]] = row[2]

## jgi_summarize_bam_contig_depths output (binning.smk::contig_abundance):
## header contigName, contigLen, totalAvgDepth, <bam>, <bam>-var. contig_id -> cov
coverage_by_contig = {}
with open(snakemake.input.coverage) as fh:
    reader = csv.reader(fh, delimiter="\t")
    next(reader)  # header
    for row in reader:
        if len(row) >= 3:
            coverage_by_contig[row[0]] = row[2]

## barrnap GFF (binqc.smk::barrnap_predict, bac+arc union): bin_id -> {5S, 16S, 23S found}
rrna_types = {}
for rrna_path in snakemake.input.rrna:
    bin_id = os.path.basename(rrna_path).removesuffix(".rrna.gff")
    with open(rrna_path) as fh:
        content = fh.read()
    rrna_types[bin_id] = {
        "5S": "5S_rRNA" in content,
        "16S": "16S_rRNA" in content,
        "23S": "23S_rRNA" in content,
    }

## tRNAscan-SE output (binqc.smk::trnascan_predict): bin_id -> tRNA gene count.
## Standard legacy tabular output has a 3-line header before data rows.
trna_count = {}
for trna_path in snakemake.input.trna:
    bin_id = os.path.basename(trna_path).removesuffix(".trna.txt")
    with open(trna_path) as fh:
        lines = fh.readlines()
    trna_count[bin_id] = sum(1 for line in lines[3:] if line.strip())

## GUNC (binqc.smk::gunc_run): bin_id -> clade separation score.
gunc_css = {}
matches = glob.glob(os.path.join(snakemake.input.gunc_report, "*maxCSS_level.tsv"))
if matches:
    with open(matches[0]) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            bin_id = row["genome"].removesuffix(".fasta")
            gunc_css[bin_id] = row.get("clade_separation_score", "NA")

min_completeness = float(snakemake.params.mimag_min_completeness)
max_contamination = float(snakemake.params.mimag_max_contamination)
min_trna = int(snakemake.params.mimag_min_trna)


def mimag_hq(bin_id, completeness, contamination):
    try:
        completeness, contamination = float(completeness), float(contamination)
    except ValueError:
        return "NA"
    types = rrna_types.get(bin_id, {})
    has_rrna = types.get("5S") and types.get("16S") and types.get("23S")
    has_trna = trna_count.get(bin_id, 0) >= min_trna
    is_hq = completeness > min_completeness and contamination < max_contamination and has_rrna and has_trna
    return "Y" if is_hq else "N"


summary_rows = []
for bin_id in passed_bins:
    src_fasta = os.path.join(snakemake.input.passed_bins, f"{bin_id}.fasta")
    dst_fasta = os.path.join(snakemake.output.finalbins, f"{bin_id}.fasta")

    with open(src_fasta) as src, open(dst_fasta, "w") as dst:
        for line in src:
            if line.startswith(">"):
                contig_id = line[1:].split()[0].strip()
                cov = coverage_by_contig.get(contig_id, "NA")
                circ = circular_by_contig.get(contig_id, "NA")
                dst.write(f">{contig_id} cov={cov} circular={circ}\n")
            else:
                dst.write(line)

    bin_stats = checkm2.get(bin_id, {})
    completeness = bin_stats.get("completeness", "NA")
    contamination = bin_stats.get("contamination", "NA")
    genome_size = bin_stats.get("genome_size", "NA")
    n_contig = bin_stats.get("n_contig", "NA")
    n50 = bin_stats.get("n50", "NA")
    gc = bin_stats.get("gc", "NA")
    summary_rows.append(
        [
            snakemake.wildcards.sample,
            snakemake.wildcards.assembler,
            bin_id,
            completeness,
            contamination,
            "NA",
            genome_size,
            n_contig,
            n50,
            gc,
            taxa.get(bin_id, "NA"),
            gunc_css.get(bin_id, "NA"),
            mimag_hq(bin_id, completeness, contamination),
        ]
    )

with open(snakemake.output.summary, "w") as fh:
    fh.write(
        "#Sample\tAssembler\tBinID\tCompleteness\tContamination\tStrain heterogeneity\t"
        "GenomeSize(bp)\tN_Contig\tN50(bp)\tGC\tGTDB-Taxa\tGUNC_CSS\tMIMAG_HQ\n"
    )
    for row in summary_rows:
        fh.write("\t".join(row) + "\n")
