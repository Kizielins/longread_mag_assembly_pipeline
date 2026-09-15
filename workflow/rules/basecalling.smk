"""Basecalling/: optional Dorado basecalling stage, upstream of everything
else in this pipeline. Off by default
Only runs for samples whose samples_sheet row sets `pod5_dir` instead of
`fastq`
"""

def raw_reads_path(wildcards):
    """Per-sample raw reads feeding 00-QC (nanoplot_raw, chopper_filter) -
    either the user-supplied fastq directly, or this stage's basecalled
    output if that sample's samples_sheet row set `pod5_dir` instead."""
    if wildcards.sample in SAMPLES_POD5:
        return f"{OUTDIR}/Basecalling/{wildcards.sample}/basecalled.fastq.gz"
    return SAMPLES_FASTQ[wildcards.sample]


def _dorado_model_string():
    mods = config.get("dorado_modified_bases", "").strip()
    base = config["dorado_model"]
    return f"{base},{mods}" if mods else base


def _dorado_pod5_input(wildcards):
    """pod5_dir is expected to be the parent directory containing pod5_pass/
    and (optionally) pod5_fail/ subdirectories recursive picks up both when 
    dorado_include_fail is set"""
    base = SAMPLES_POD5[wildcards.sample]
    if cfg_bool("dorado_include_fail", True):
        return base
    return os.path.join(base, "pod5_pass")


rule dorado_basecall:
    input:
        pod5_dir=_dorado_pod5_input,
    output:
        bam=f"{OUTDIR}/Basecalling/{{sample}}/basecalled.bam",
        fastq=f"{OUTDIR}/Basecalling/{{sample}}/basecalled.fastq.gz",
    params:
        bin=expand_path(config.get("dorado_bin", "dorado")),
        model=_dorado_model_string(),
        device=config.get("dorado_device", "cuda:all"),
    threads: config["threads"]
    resources:
        gpu=1,
    log:
        f"{OUTDIR}/Basecalling/{{sample}}/dorado.log",
    shell:
        """
        {params.bin} basecaller {params.model} {input.pod5_dir} --recursive \
            --device {params.device} > {output.bam} 2> {log}
        samtools fastq {output.bam} 2>> {log} | gzip > {output.fastq}
        """
