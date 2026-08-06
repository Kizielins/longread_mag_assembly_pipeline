"""06-StrainPhasing: Track B, strain-level phasing via Strainy (Kazantseva
et al., Nat Methods 2024), plus a second, independent source of
strain-resolved genomes when myloasm is one of Track A's assemblers.

Inputs:
- assembly graph - whichever Track A assembler's own assembly_graph.gfa is
  picked by _strainy_assembler() below 
- reads - QC'd reads from qc.smk

Outputs:
- strainy-result/ - Strainy's phased-strain output
- myloasm-alternate-contigs.fa - myloasm's own near-duplicate/strain-level
  contigs, excluded from assembly_primary.fa 
"""


def _strainy_ref_flag(wildcards, input):
    """metaFlye: real graph, Strainy's tested input shape. nanoMDBG/myloasm:
    fasta_ref instead of their own graph"""
    if wildcards.assembler == "metaflye":
        return f"--gfa_ref {input.gfa}"
    return f"--fasta_ref {input.fasta}"


if cfg_bool("run_strain_phasing"):

    rule strain_phasing:
        input:
            gfa=f"{OUTDIR}/01-LongAssemblies/{{sample}}/{{assembler}}/assembly_graph.gfa",
            fasta=f"{OUTDIR}/01-LongAssemblies/{{sample}}/{{assembler}}/assembly.fasta",
            reads=qc_reads_path,
        output:
            directory(f"{OUTDIR}/06-StrainPhasing/{{sample}}/{{assembler}}/strainy-result"),
        params:
            ref_flag=_strainy_ref_flag,
        threads: config["threads"]
        conda:
            config["main_env"]
        log:
            f"{OUTDIR}/06-StrainPhasing/{{sample}}/{{assembler}}/strainy.log",
        shell:
            """
            strainy {params.ref_flag} --fastq {input.reads} --mode nano \
                --output {output} --threads {threads} > {log} 2>&1
            """


if "myloasm" in ASSEMBLERS:

    rule myloasm_alternate_contigs:
        """Publishes myloasm's own near-duplicate/strain-level contigs
        (assembly.smk::assemble_myloasm's `alternate` output) under
        06-StrainPhasing/ as a second, independent source of strain-resolved
        genome"""
        input:
            rules.assemble_myloasm.output.alternate,
        output:
            f"{OUTDIR}/06-StrainPhasing/{{sample}}/myloasm/alternate-contigs.fa",
        shell:
            "cp {input} {output}"
