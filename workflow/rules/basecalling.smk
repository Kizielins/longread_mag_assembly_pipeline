"""Basecalling/: optional Dorado basecalling stage, upstream of everything
else in this pipeline. Off by default
Only runs for samples whose samples_sheet row sets `pod5_dir` instead of
`fastq`
"""

def raw_reads_path(wildcards):
    """Per-sample raw reads feeding 00-QC (nanoplot_raw, chopper_filter) -
    either the user-supplied fastq directly, or this stage's basecalled
    output if that sample's samples_sheet row set `pod5_dir` instead."""
    if wildcards.sample in SAMPLES_POD5_BARCODED:
        run_id, barcode = SAMPLES_POD5_BARCODED[wildcards.sample]
        return f"{OUTDIR}/Basecalling/{run_id}/demux/{barcode}.fastq.gz"
    return SAMPLES_FASTQ[wildcards.sample]


def _dorado_modbase_flag():
    mods = config.get("dorado_modified_bases", "").strip()
    return f"--modified-bases-models {mods}" if mods else ""

def _dorado_pod5_input(pod5_dir):
    """pod5_dir is expected to be the parent directory containing pod5_pass/
    and (optionally) pod5_fail/ subdirectories recursive picks up both when 
    dorado_include_fail is set"""
    if cfg_bool("dorado_include_fail", True):
        return pod5_dir
    return os.path.join(pod5_dir, "pod5_pass")


rule dorado_basecall_run:
    input:
        pod5_dir=lambda wc: _dorado_pod5_input(RUNS[wc.run]),
    output:
        bam=f"{OUTDIR}/Basecalling/{{run}}/basecalled.bam",
    params:
        bin=expand_path(config.get("dorado_bin", "dorado")),
        model=config["dorado_model"],
        modbase_flag=_dorado_modbase_flag(),
        device=config.get("dorado_device", "cuda:all"),
        kit=config["dorado_barcode_kit"],
    threads: config["threads"]
    resources:
        gpu=1,
    log:
        f"{OUTDIR}/Basecalling/{{run}}/dorado.log",
    shell:
        """
        {params.bin} basecaller {params.model} {input.pod5_dir} --recursive \
            --kit-name {params.kit} {params.modbase_flag} --trim none --device {params.device} > {output.bam} 2> {log}
        """

rule dorado_demux:
    input:
        bam=rules.dorado_basecall_run.output.bam,
    output:
        outdir=directory(f"{OUTDIR}/Basecalling/{{run}}/demux-raw"),
    params:
        bin=expand_path(config.get("dorado_bin", "dorado")),
    threads: config["threads"]
    log:
        f"{OUTDIR}/Basecalling/{{run}}/demux.log",
    shell:
        """
        {params.bin} demux --no-classify --output-dir {output.outdir} {input.bam} > {log} 2>&1
        """

rule barcode_fastq:
    input:
        demux_dir=rules.dorado_demux.output.outdir,
    output:
        f"{OUTDIR}/Basecalling/{{run}}/demux/{{barcode}}.fastq.gz",
    log:
        f"{OUTDIR}/Basecalling/{{run}}/demux/{{barcode}}.fastq.log",
    shell:
        """
        bc_bam=$(find {input.demux_dir} -iname "*{wildcards.barcode}*.bam" | head -n1)
        if [ -z "$bc_bam" ]; then
            echo "No demuxed BAM found for {wildcards.barcode} under {input.demux_dir}" > {log}
            exit 1
        fi
        samtools fastq "$bc_bam" 2> {log} | gzip > {output}
        """
