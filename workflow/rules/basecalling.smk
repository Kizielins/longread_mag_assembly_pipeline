"""Basecalling/: optional Dorado basecalling stage, upstream of everything
else in this pipeline. Off by default
Only runs for samples whose samples_sheet row sets `pod5_dir` instead of
`fastq`
"""


wildcard_constraints:
    barcode=r"barcode\d{2}",


def raw_reads_path(wildcards):
    """Per-sample raw reads feeding 00-QC (nanoplot_raw, chopper_filter)."""
    if wildcards.sample in SAMPLES_POD5_BARCODED:
        run_id, barcode = SAMPLES_POD5_BARCODED[wildcards.sample]
        return f"{OUTDIR}/Basecalling/{run_id}/demux/{barcode}.fastq.gz"
    return SAMPLES_FASTQ[wildcards.sample]


def modbam_path(wildcards):
    """Per-sample BAM WITH modified-base tags, for methylation analysis.

    Deliberately separate from raw_reads_path: the fastq branch loses
    MM/ML. Returns None for plain-fastq samples"""
    if wildcards.sample in SAMPLES_POD5_BARCODED:
        run_id, barcode = SAMPLES_POD5_BARCODED[wildcards.sample]
        return f"{OUTDIR}/Basecalling/{run_id}/demux/{barcode}.bam"
    return None


def _dorado_modbase_flag():
    """dorado_modified_bases holds full modbase model names (comma-
    separated, e.g. "dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1,...")
    """
    mods = config.get("dorado_modified_bases", "").strip()
    return f"--modified-bases-models {mods}" if mods else ""


def _dorado_pod5_input(pod5_dir):
    """pod5_dir is the run's parent directory, holding only pod5_pass/ and
    (when present) pod5_fail/ -- no pod5_skip/ or other sibling pod5 tree
    in this pipeline's real data, so --recursive is safe to point at the
    parent as a whole rather than staging pod5_pass/pod5_fail explicitly.
    If a pod5_skip/ (or anything else pod5-bearing) can ever show up
    alongside them, --recursive would sweep that in too -- revisit then.
    """
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
            --kit-name {params.kit} {params.modbase_flag} \
            --device {params.device} > {output.bam} 2> {log}
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
        {params.bin} demux --no-classify --emit-summary \
            --threads {threads} \
            --output-dir {output.outdir} {input.bam} > {log} 2>&1
        """

rule barcode_bam:
    input:
        demux_dir=rules.dorado_demux.output.outdir,
    output:
        bam=f"{OUTDIR}/Basecalling/{{run}}/demux/{{barcode}}.bam",
    log:
        f"{OUTDIR}/Basecalling/{{run}}/demux/{{barcode}}.bam.log",
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})"

        mapfile -t bc_dirs < <(find {input.demux_dir} -type d -name "{wildcards.barcode}" | sort)
        if [ "${{#bc_dirs[@]}}" -eq 0 ]; then
            echo "dorado demux classified no reads to {wildcards.barcode} in run {wildcards.run}" \
                 | tee {log} >&2
            exit 1
        fi
        if [ "${{#bc_dirs[@]}}" -gt 1 ]; then
            echo "Multiple demux dirs matched {wildcards.barcode} in run {wildcards.run} (e.g. bam_pass/ + bam_fail/?) -- not yet handled, extend this rule once seen for real:" \
                 | tee {log} >&2
            printf '%s\n' "${{bc_dirs[@]}}" | tee -a {log} >&2
            exit 1
        fi

        mapfile -t bc_bams < <(find "${{bc_dirs[0]}}" -maxdepth 1 -type f -name '*.bam' | sort)
        if [ "${{#bc_bams[@]}}" -ne 1 ]; then
            echo "Expected exactly one BAM in ${{bc_dirs[0]}} for {wildcards.barcode} in run {wildcards.run}, found ${{#bc_bams[@]}} -- not yet handled, extend this rule once seen for real:" \
                 | tee {log} >&2
            printf '%s\n' "${{bc_bams[@]}}" | tee -a {log} >&2
            exit 1
        fi

        echo "Linking ${{bc_bams[0]}} -> {output.bam}" > {log}
        ln -f "${{bc_bams[0]}}" {output.bam} 2>/dev/null || cp "${{bc_bams[0]}}" {output.bam}
        """


rule barcode_fastq:
    input:
        bam=rules.barcode_bam.output.bam,
    output:
        f"{OUTDIR}/Basecalling/{{run}}/demux/{{barcode}}.fastq.gz",
    params:
        gzip=config.get("gzip_bin", "pigz"),
    threads: 4
    log:
        f"{OUTDIR}/Basecalling/{{run}}/demux/{{barcode}}.fastq.log",
    shell:
        """
        samtools fastq -@ {threads} {input.bam} 2> {log} \
          | {params.gzip} -p {threads} > {output}
        """


