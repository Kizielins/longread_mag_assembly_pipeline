"""03-LongBins/BIN_REFINEMENT: DAS_Tool bin refinement. 

Inputs:
- per-binner bin fastas - metabat2/SemiBin2(/COMEBin) output from binning.smk
- assembly - polished assembly

Outputs:
- scaffolds2bin/{binner}.tsv - contig->bin table per binner, DAS_Tool's expected input format
- refined-bins/ - DAS_Tool's consensus bin set, renamed to bin.1.fa, bin.2.fa, ... 
"""


def _scaffolds2bin_input(wc):
    mapping = {
        "metabat2": rules.metabat2_binning.output.bins,
        "semibin2": rules.semibin2_binning.output.bins,
    }
    if cfg_bool("use_comebin"):
        mapping["comebin"] = rules.comebin_binning.output.bins
    return mapping[wc.binner]


rule scaffolds2bin:
    input:
        bins=_scaffolds2bin_input,
    output:
        f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/BIN_REFINEMENT/dastool/scaffolds2bin/{{binner}}.tsv",
    params:
        ext="fa",
    conda:
        config["main_env"]
    shell:
        "Fasta_to_Contig2Bin.sh -i {input.bins} -e {params.ext} > {output}"


checkpoint das_tool_refine:
    input:
        tsvs=lambda wc: expand(
            f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/BIN_REFINEMENT/dastool/scaffolds2bin/{{binner}}.tsv",
            sample=wc.sample,
            assembler=wc.assembler,
            binner=active_binners(),
        ),
        assembly=rules.medaka_polish_assembly.output.consensus,
    output:
        dastool_bins=directory(f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/BIN_REFINEMENT/dastool/dastool_DASTool_bins"),
        bins_dir=directory(f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/BIN_REFINEMENT/refined-bins"),
    params:
        outprefix=f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/BIN_REFINEMENT/dastool/dastool",
        score_threshold=config["dastool_score_threshold"],
        labels=",".join(active_binners()),
        tsv_list=lambda wc, input: ",".join(input.tsvs),
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/BIN_REFINEMENT/bin_refinement.log",
    shell:
        """
        DAS_Tool -i {params.tsv_list} -l {params.labels} -c {input.assembly} \
            -o {params.outprefix} --write_bins --score_threshold {params.score_threshold} \
            --threads {threads} > {log} 2>&1

        mkdir -p {output.bins_dir}
        i=1
        for f in {output.dastool_bins}/*.fa; do
            [ -f "$f" ] || continue
            cp "$f" "{output.bins_dir}/bin.${{i}}.fa"
            i=$((i + 1))
        done
        ls {output.bins_dir}/bin.*.fa >/dev/null 2>&1
        """
