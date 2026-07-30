"""03-LongBins/INITIAL_BINNING: read mapping -> coverage, then metabat2 +
SemiBin2 (+ optional COMEBin) binning.
"""


rule map_reads_to_assembly:
"""Sorted, indexed BAM of long reads mapped back to the polished assembly."""
    input:
        assembly=rules.medaka_polish_assembly.output.consensus,
        reads=qc_reads_path,
    output:
        bam=f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/mapping/contig.mapped.sorted.bam",
        bai=f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/mapping/contig.mapped.sorted.bam.bai",
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/mapping/mapping.log",
    shell:
        """
        asm_gb=$(awk '/^>/{{next}}{{n+=length($0)}}END{{print int(n/1000000000)+1}}' {input.assembly}) # assembly size in GB
        mini_ref=$(( asm_gb > 4 ? asm_gb : 4 ))
        minimap2 -a -x map-ont -I${{mini_ref}}g -t {threads} {input.assembly} {input.reads} 2> {log} \
            | samtools view -h -b -S -@ {threads} - \
            | samtools view -b -F 4 -@ {threads} - \ 
            | samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """


rule contig_abundance:
    """Per-contig depth via jgi_summarize_bam_contig_depths"""
    input:
        bam=rules.map_reads_to_assembly.output.bam,
        bai=rules.map_reads_to_assembly.output.bai,
    output:
        f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/depth.txt",
    conda:
        config["main_env"]
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam}"


rule metabat2_binning:
    input:
        assembly=rules.medaka_polish_assembly.output.consensus,
        depth=rules.contig_abundance.output[0],
    output:
        bins=directory(f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/metabat2/metabat2-bins"),
        log=f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/metabat2/bin.log",
    threads: config["threads"]
    conda:
        config["main_env"]
    shell:
        """
        metabat2 -t {threads} -i {input.assembly} -o {output.bins}/bin -a {input.depth} > {output.log}
        grep -q "formed" {output.log}
        """


rule semibin2_binning:
    """ Uses --depth-metabat2 to reuse contig_abundance's depth.txt directly
    instead of -b/BAM """ 
    input:
        assembly=rules.medaka_polish_assembly.output.consensus,
        depth=rules.contig_abundance.output[0],
    output:
        bins=directory(f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/semibin/semibin-output/output_bins"),
    params:
        environment=config["semibin_env"],
        outdir=f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/semibin/semibin-output",
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/semibin/semibin.log",
    shell:
        """
        SemiBin2 single_easy_bin --environment {params.environment} -i {input.assembly} --depth-metabat2 {input.depth} \
            -p {threads} --sequencing-type=long_read --compression none -o {params.outdir} > {log} 2>&1
        ls {output.bins}/*.fa >/dev/null 2>&1
        """


if cfg_bool("use_comebin"):

    rule comebin_binning:
        """Optional third binner for extra ensemble diversity into DAS_Tool.
        Needs its own env (config["comebin_env"] """

        input:
            assembly=rules.medaka_polish_assembly.output.consensus,
            bam=rules.map_reads_to_assembly.output.bam,
            bai=rules.map_reads_to_assembly.output.bai,
        output:
            bins=directory(f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/comebin/comebin-bins"),
        params:
            bamdir=f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/comebin/bam",
            outdir=f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/comebin/comebin-output",
        threads: config["threads"]
        conda:
            config["comebin_env"]
        log:
            f"{OUTDIR}/03-LongBins/{{sample}}/{{assembler}}/INITIAL_BINNING/comebin/comebin.log",
        shell:
            """
            mkdir -p {params.bamdir}
            ln -sf $(realpath {input.bam}) {params.bamdir}/
            run_comebin.sh -a {input.assembly} -p {params.bamdir} -o {params.outdir} -t {threads} > {log} 2>&1
            mkdir -p {output.bins}
            cp {params.outdir}/comebin_res/comebin_res_bins/*.fa {output.bins}/ 2>/dev/null \
                || find {params.outdir}/comebin_res -name '*.fa' -exec cp {{}} {output.bins}/ \; 2>/dev/null || true
            ls {output.bins}/*.fa >/dev/null 2>&1
            """
