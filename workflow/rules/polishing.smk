"""02-Polishing: single-round medaka polishing of the whole assembly.

Inputs:
- draft - raw assembly from assembly.smk
- reads - QC'd reads from qc.smk

Outputs:
- polished_assembly.fasta -- medaka consensus assembly
"""

rule medaka_polish_assembly:
    input:
        draft=assembly_output,
        reads=qc_reads_path,
    output:
        consensus=f"{OUTDIR}/02-Polishing/{{sample}}/{{assembler}}/polished_assembly.fasta",
    params:
        outdir=f"{OUTDIR}/02-Polishing/{{sample}}/{{assembler}}/medaka",
        model=config["medaka_model"],
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/02-Polishing/{{sample}}/{{assembler}}/medaka.log",
    shell:
        """
        medaka_consensus -m {params.model} -i {input.reads} -d {input.draft} -o {params.outdir} -t {threads} > {log} 2>&1
        cp {params.outdir}/consensus.fasta {output.consensus}
        rm -rf {params.outdir}/calls_to* {params.outdir}/consensus_probs.hdf {params.outdir}/consensus.fasta.gaps_in_draft_coords.bed
        """
