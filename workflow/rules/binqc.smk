"""Bin QC: CheckM2 completeness/contamination, GUNC chimera detection
 and quality filtering (>=50% complete, <=10% contaminated) down to the final MAG set.
Filtering is two stages: filter_checkm2_bins (completeness/contamination)
feeds gunc_run, which only scans survivors, which feeds filter_quality_bins.

Inputs:
- refined bins - refinement.smk's refined-bins/{bin}.fa, one per bin ID from get_bin_ids()

Outputs:
- polished-bins/{bin}.fasta - staged per-bin fasta copy, common input for the QC tools below
- checkm2-result/quality_report.tsv - completeness/contamination + genome size/N50/GC/contig count per bin
- checkm2-passed-bins/, checkm2_passed_bins.list - stage-1 survivors
- gunc-result/ - chimerism check on stage-1 survivors
- {bin}.rrna.gff / {bin}.trna.txt - barrnap/tRNAscan-SE predictions, full refined set
- passed_bins.list / passed-bins/ - final survivors, feeding GTDB-Tk
"""

rule prepare_bin_fasta:
    """ Just for input clarity"""
    input:
        f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/BIN_REFINEMENT/refined-bins/{{bin}}.fa",
    output:
        f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.fasta",
    shell:
        "ln -sf $(realpath {input}) {output}"


rule checkm2_predict:
    input:
        bins=lambda wc: expand(
            f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.fasta",
            sample=wc.sample,
            assembler=wc.assembler,
            bin=get_bin_ids(wc),
        ),
    output:
        f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/checkm2-result/quality_report.tsv",
    params:
        indir=f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins",
        outdir=f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/checkm2-result",
        db=expand_path(config["checkm2_db"]),
    threads: config["threads"]
    conda:
        config["checkm2_env"]
    log:
        f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/checkm2-result.log",
    shell:
        "checkm2 predict --input {params.indir} --output-directory {params.outdir} "
        "--database_path {params.db} --extension fasta --threads {threads} --force > {log} 2>&1"

rule filter_checkm2_bins:
    """Completeness/contamination filter"""
    input:
        checkm2_report=rules.checkm2_predict.output[0],
        bins=lambda wc: expand(
            f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.fasta",
            sample=wc.sample,
            assembler=wc.assembler,
            bin=get_bin_ids(wc),
        ),
    output:
        passed_list=f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/checkm2_passed_bins.list",
        finalbins=directory(f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/checkm2-passed-bins"),
    params:
        min_completeness=config["min_completeness"],
        max_contamination=config["max_contamination"],
    script:
        "../scripts/filter_checkm2_bins.py"

rule gunc_run:
    input:
        bins=rules.filter_checkm2_bins.output.finalbins,
    output:
        directory(f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/gunc-result"),
    params:
        db=expand_path(config["gunc_db"]),
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/gunc-result.log",
    shell:
        "mkdir -p {output} && "
        "gunc run --input_dir {input.bins} --file_suffix .fasta --out_dir {output} --db_file {params.db} "
        "--threads {threads} > {log} 2>&1"

rule barrnap_predict:
    """5S/16S/23S rRNA detection, for MIMAG high-quality status. """
    input:
        rules.prepare_bin_fasta.output[0],
    output:
        f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.rrna.gff",
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.barrnap.log",
    shell:
        """
        {{ barrnap --kingdom bac --quiet {input} || true; barrnap --kingdom arc --quiet {input} || true; }} \
            > {output} 2> {log}
        """


rule trnascan_predict:
    """tRNA gene count, for MIMAG high-quality status (>=18 tRNAs)."""
    input:
        rules.prepare_bin_fasta.output[0],
    output:
        f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.trna.txt",
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.trnascan.log",
    shell:
        "tRNAscan-SE -G --forceow -o {output} -m /dev/null {input} > {log} 2>&1"


rule filter_quality_bins:
    """Drops GUNC-flagged chimeras from filter_checkm2_bins' survivors, and
    copies the final FASTAs into staging/passed-bins/ for GTDB-Tk."""
    input:
        checkm2_passed_list=rules.filter_checkm2_bins.output.passed_list,
        bins=rules.filter_checkm2_bins.output.finalbins,
        gunc_report=rules.gunc_run.output[0],
    output:
        passed_list=f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/passed_bins.list",
        finalbins=directory(f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/passed-bins"),
    script:
        "../scripts/filter_quality_bins.py"
