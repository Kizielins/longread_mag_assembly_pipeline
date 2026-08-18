"""GTDB-Tk classification of the quality-passed bins (05-GTDBTK)"""

rule gtdbtk_classify:
    """GTDB-Tk has no CLI flag for its reference data path -- it strictly
    requires the GTDBTK_DATA_PATH environment variable. Set inline in this
    rule's own shell command (from config["gtdbtk_data_path"]) rather than
    relying on the calling shell having it exported, so this rule is
    self-contained regardless of how snakemake itself was invoked."""
    input:
        genome_dir=rules.filter_quality_bins.output.finalbins,
        passed_list=rules.filter_quality_bins.output.passed_list,
    output:
        result=directory(f"{OUTDIR}/05-GTDBTK/{{sample}}/{{assembler}}/gtdbtk-result"),
        taxa=f"{OUTDIR}/05-GTDBTK/{{sample}}/{{assembler}}/tmp.taxa",
    params:
        data_path=expand_path(config["gtdbtk_data_path"]),
        pplacer_cpus=config["gtdbtk_pplacer_cpus"],
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/05-GTDBTK/{{sample}}/{{assembler}}/gtdbtk-result.log",
    shell:
        """
        GTDBTK_DATA_PATH={params.data_path} gtdbtk classify_wf --genome_dir {input.genome_dir} -x fasta \
            --out_dir {output.result} --cpus {threads} --pplacer_cpus {params.pplacer_cpus} --skip_ani_screen > {log} 2>&1
        cat {output.result}/classify/gtdbtk.*summary.tsv | grep -v '^user_genome' \
            | awk -F"\\t" '{{print $1"\\t"$2}}' > {output.taxa}
        """


rule build_genome_summary:
    input:
        checkm2_report=rules.checkm2_predict.output[0],
        rrna=lambda wc: expand(
            f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.rrna.gff",
            sample=wc.sample,
            assembler=wc.assembler,
            bin=get_bin_ids(wc),
        ),
        trna=lambda wc: expand(
            f"{OUTDIR}/04-BinQC/{{sample}}/{{assembler}}/staging/polished-bins/{{bin}}.trna.txt",
            sample=wc.sample,
            assembler=wc.assembler,
            bin=get_bin_ids(wc),
        ),
        passed_bins=rules.filter_quality_bins.output.finalbins,
        passed_list=rules.filter_quality_bins.output.passed_list,
        gunc_report=rules.gunc_run.output[0],
        taxa=rules.gtdbtk_classify.output.taxa,
        contigs_info=contigs_info_output,
        coverage=rules.contig_abundance.output[0],
    output:
        summary=f"{OUTDIR}/{{sample}}.{{assembler}}.pipeline.ont.genome.summary",
        finalbins=directory(f"{OUTDIR}/Final-bins/{{sample}}/{{assembler}}"),
    params:
        mimag_min_completeness=config["mimag_min_completeness"],
        mimag_max_contamination=config["mimag_max_contamination"],
        mimag_min_trna=config["mimag_min_trna"],
    script:
        "../scripts/build_genome_summary.py"
