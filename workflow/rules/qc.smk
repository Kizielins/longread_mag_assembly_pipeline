"""00-QC: NanoPlot QC visualization, Chopper length/quality filtering,
host-read removal with hostile.

Inputs:
- raw reads

Outputs:
- nanoplot_raw/ - NanoPlot QC report on the raw, unfiltered reads
- filtered.reads.fq.gz - Chopper length/quality-filtered reads (temp, intermediate)
- dehosted.reads.fq.gz - filtered reads with human host reads removed (skipped if skip_host_removal)
- nanoplot_filtered/ - NanoPlot QC report on the final filtered reads
"""

def qc_reads_path(wildcards):
    """Final QC'd read file that assembly uses"""
    if cfg_bool("skip_host_removal"):
        return f"{OUTDIR}/00-QC/{wildcards.sample}/filtered.reads.fq.gz"
    return f"{OUTDIR}/00-QC/{wildcards.sample}/dehosted.reads.fq.gz"


rule nanoplot_raw:
    input:
        lambda wc: SAMPLES[wc.sample],
    output:
        directory(f"{OUTDIR}/00-QC/{{sample}}/nanoplot_raw"),
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/00-QC/{{sample}}/nanoplot_raw.log",
    shell:
        "NanoPlot --fastq {input} --outdir {output} --threads {threads} "
        "--loglength --N50 > {log} 2>&1"


rule chopper_filter:
    input:
        lambda wc: SAMPLES[wc.sample],
    output:
        temp(f"{OUTDIR}/00-QC/{{sample}}/filtered.reads.fq.gz"),
    params:
        min_q=config["chopper_min_quality"],
        min_l=config["chopper_min_length"],
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/00-QC/{{sample}}/chopper.log",
    shell:
        "zcat -f {input} | chopper -q {params.min_q} -l {params.min_l} --threads {threads} "
        "2> {log} | gzip > {output}"


rule remove_host_reads:
    input:
        f"{OUTDIR}/00-QC/{{sample}}/filtered.reads.fq.gz",
    output:
        f"{OUTDIR}/00-QC/{{sample}}/dehosted.reads.fq.gz",
    params:
        index_flag=lambda wc: f"--index {config['host_index']}" if config.get("host_index") else "",
        tmpdir=f"{OUTDIR}/00-QC/{{sample}}/.hostile_tmp",
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/00-QC/{{sample}}/hostile.log",
    shell:
        """
        mkdir -p {params.tmpdir}
        hostile clean --fastq1 {input} --aligner minimap2 {params.index_flag} \
            --output {params.tmpdir} --threads {threads} > {log} 2>&1
        mv {params.tmpdir}/*.fastq.gz {output}
        rm -rf {params.tmpdir}
        """


rule nanoplot_filtered:
    input:
        qc_reads_path,
    output:
        directory(f"{OUTDIR}/00-QC/{{sample}}/nanoplot_filtered"),
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/00-QC/{{sample}}/nanoplot_filtered.log",
    shell:
        "NanoPlot --fastq {input} --outdir {output} --threads {threads} "
        "--loglength --N50 > {log} 2>&1"
