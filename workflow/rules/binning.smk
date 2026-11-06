"""03-LongBins/INITIAL_BINNING: read mapping -> coverage, then metabat2 +
SemiBin2 (+ optional COMEBin) binning.

Inputs:
- assembly - polished assembly from polishing.smk
- reads - QC'd reads from qc.smk

Outputs:
- contig.mapped.sorted.bam(.bai) - reads mapped back onto the assembly
- depth.txt - per-contig coverage (jgi_summarize_bam_contig_depths)
- metabat2-bins/ - metabat2's bin fastas
- output_bins/ - SemiBin2's bin fastas
- comebin-bins/ - COMEBin's bin fastas (only if use_comebin is set)

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
        "jgi_summarize_bam_contig_depths --percentIdentity 70 --outputDepth {output} {input.bam}"


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
        mkdir -p {output.bins}
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
        mkdir -p {output.bins}
        if ! ls {output.bins}/*.fa >/dev/null 2>&1; then
            echo "WARNING: SemiBin2 produced zero bins for this (sample, assembler) - contributing nothing to the DAS_Tool ensemble. See semibin.log for details." >> {log}
        fi
        """

if cfg_bool("use_comebin"):

    rule comebin_binning:
        """Optional third binner for extra ensemble diversity into DAS_Tool.
        Needs its own env (config["comebin_env"] 
        
        COMEBin can legitimately crash mid-run on some assemblies: its
        internal marker-gene seed-selection step (test_getmarker_*.pl) can
        find zero qualifying seed contigs for a given k-means ensemble
        method (observed on a real run - the script itself exits 0 and
        writes nothing, which is normal behavior for a low-marker-density
        input), but COMEBin's Python side doesn't skip that method when
        that happens - it proceeds into its own bundled CheckM v1 scoring
        step anyway, which then hard-crashes with FileNotFoundError trying
        to list a bins/ directory that was never written (upstream bug,
        unfixed as of ziyewang/COMEBin#53, same symptom reported there).
        Check comebin.log if comebin-bins/ turns up empty"""

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
            r"""
            mkdir -p {params.bamdir} {output.bins}
            ln -sf $(realpath {input.bam}) {params.bamdir}/
            run_comebin.sh -a {input.assembly} -p {params.bamdir} -o {params.outdir} -t {threads} > {log} 2>&1 \
                || echo "WARNING: run_comebin.sh exited non-zero -- COMEBin is optional/best-effort, continuing with whatever bins (if any) it produced instead of failing this (sample, assembler)'s DAS_Tool refinement (see docstring)." >> {log}
            cp {params.outdir}/comebin_res/comebin_res_bins/*.fa {output.bins}/ 2>/dev/null \
                || find {params.outdir}/comebin_res -name '*.fa' -exec cp {{}} {output.bins}/ \; 2>/dev/null || true
            if ! ls {output.bins}/*.fa >/dev/null 2>&1; then
                echo "WARNING: COMEBin produced zero usable bins for this (sample, assembler) -- contributing nothing to the DAS_Tool ensemble. See comebin.log for details." >> {log}
            fi
            """
