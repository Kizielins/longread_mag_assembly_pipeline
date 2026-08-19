"""01-LongAssemblies: long-read metagenome assembly.
Inputs:
- reads - QC'd reads from qc.smk

Output (regardless of assembler):
    - assembly.fasta
    - assembly_graph.gfa
    - contigs_info.tsv (seq_name, length, circular)

assembly.fasta/ assembly_graph.gfa here are the raw results, saved for further use with Track B (strain phasing).
"""

def assembly_output(wildcards):
    return f"{OUTDIR}/01-LongAssemblies/{wildcards.sample}/{wildcards.assembler}/assembly.fasta"


def contigs_info_output(wildcards):
    return f"{OUTDIR}/01-LongAssemblies/{wildcards.sample}/{wildcards.assembler}/contigs_info.tsv"

def flye_extra_flags():
    """Optional flye flags built from config/config.yaml's flye_read_error /
    flye_asm_coverage / flye_genome_size - all empty by default """
    flags = []
    read_error = str(config.get("flye_read_error", "")).strip()
    if read_error:
        flags.append(f"--read-error {read_error}")
    asm_coverage = str(config.get("flye_asm_coverage", "")).strip()
    genome_size = str(config.get("flye_genome_size", "")).strip()
    if asm_coverage and genome_size:
        flags.append(f"--asm-coverage {asm_coverage} -g {genome_size}")
    return " ".join(flags)

def myloasm_extra_flags():
    """Optional myloasm flags built from config/config.yaml's
    myloasm_quality_value_cutoff / myloasm_bloom_filter_size /
    myloasm_dereplication_ani / myloasm_dereplication_length """
    flags = []
    quality_cutoff = str(config.get("myloasm_quality_value_cutoff", "")).strip()
    if quality_cutoff:
        flags.append(f"--quality-value-cutoff {quality_cutoff}")
    bloom_size = str(config.get("myloasm_bloom_filter_size", "")).strip()
    if bloom_size:
        flags.append(f"--bloom-filter-size {bloom_size}")
    dereplication_ani = str(config.get("myloasm_dereplication_ani", "")).strip()
    dereplication_length = str(config.get("myloasm_dereplication_length", "")).strip()
    if dereplication_ani and dereplication_length:
        flags.append(f"--dereplication-ani {dereplication_ani} --dereplication-length {dereplication_length}")
    return " ".join(flags)


rule assemble_metaflye:
    """`-i 0 --keep-haplotypes --no-alt-contigs`: Flye's own internal
    polishing is disabled (`-i 0`) and haplotype bubbles are kept in the
    graph - required for Track B """
    input:
        reads=qc_reads_path,
    output:
        assembly=f"{OUTDIR}/01-LongAssemblies/{{sample}}/metaflye/assembly.fasta",
        contigs_info=f"{OUTDIR}/01-LongAssemblies/{{sample}}/metaflye/contigs_info.tsv",
        graph=f"{OUTDIR}/01-LongAssemblies/{{sample}}/metaflye/assembly_graph.gfa",
        log=f"{OUTDIR}/01-LongAssemblies/{{sample}}/metaflye/flye.log",
        assembly_info=f"{OUTDIR}/01-LongAssemblies/{{sample}}/metaflye/assembly_info.txt",
        graph_path=f"{OUTDIR}/01-LongAssemblies/{{sample}}/metaflye/contig_graph_path.tsv",
    params:
        tmpdir=f"{OUTDIR}/01-LongAssemblies/{{sample}}/metaflye/tmp",
        extra_flags=flye_extra_flags(),
    threads: lambda wc: min(config["threads"], 128)  # flye caps out around 128 threads
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/01-LongAssemblies/{{sample}}/metaflye/flye.log.debug",
    shell:
        """
        mkdir -p {params.tmpdir}
        flye --meta --nano-hq {input.reads} -t {threads} -i 0 --keep-haplotypes --no-alt-contigs \
            {params.extra_flags} -o {params.tmpdir} > {log} 2>&1
        grep -q "Final assembly" {params.tmpdir}/flye.log
        mv {params.tmpdir}/assembly.fasta {output.assembly}
        mv {params.tmpdir}/assembly_graph.gfa {output.graph}
        mv {params.tmpdir}/flye.log {output.log}
        mv {params.tmpdir}/assembly_info.txt {output.assembly_info}
        {{ echo -e "seq_name\\tlength\\tcircular"; tail -n+2 {output.assembly_info} \
            | awk -F"\\t" 'BEGIN{{OFS="\\t"}}{{print $1,$2,$4}}'; }} > {output.contigs_info}
        {{ echo -e "seq_name\\tgraph_path"; tail -n+2 {output.assembly_info} \
            | awk -F"\\t" 'BEGIN{{OFS="\\t"}}{{print $1,$8}}'; }} > {output.graph_path}
        rm -rf {params.tmpdir}
        """


rule assemble_nanomdbg:
    input:
        reads=qc_reads_path,
    output:
        assembly=f"{OUTDIR}/01-LongAssemblies/{{sample}}/nanomdbg/assembly.fasta",
        contigs_info=f"{OUTDIR}/01-LongAssemblies/{{sample}}/nanomdbg/contigs_info.tsv",
        graph=f"{OUTDIR}/01-LongAssemblies/{{sample}}/nanomdbg/assembly_graph.gfa",
    params:
        tmpdir=f"{OUTDIR}/01-LongAssemblies/{{sample}}/nanomdbg/tmp",
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/01-LongAssemblies/{{sample}}/nanomdbg/nanomdbg.log",
    shell:
        """
        mkdir -p {params.tmpdir}
        metaMDBG asm --out-dir {params.tmpdir} --in-ont {input.reads} --threads {threads} > {log} 2>&1
        zcat {params.tmpdir}/contigs.fasta.gz > {output.assembly}
        {{ echo -e "seq_name\\tlength\\tcircular"; grep "^>" {output.assembly} | sed 's/^>//' \
            | awk 'BEGIN{{OFS="\\t"}}{{split($2,l,"="); split($4,c,"="); \
                circ=(c[2]=="yes")?"Y":"N"; print $1,l[2],circ}}'; }} > {output.contigs_info}
        metaMDBG gfa --assembly-dir {params.tmpdir} --k 0 --threads {threads} > {params.tmpdir}/k_list.txt 2>&1
        LAST_K=$(grep -oP '(?<=- )[0-9]+' {params.tmpdir}/k_list.txt | tail -1)
        metaMDBG gfa --assembly-dir {params.tmpdir} --k $LAST_K --threads {threads} >> {log} 2>&1
        mv {params.tmpdir}/assemblyGraph_k$LAST_K.gfa {output.graph}
        rm -rf {params.tmpdir}
        """

rule assemble_myloasm:
    input:
        reads=qc_reads_path,
    output:
        assembly=f"{OUTDIR}/01-LongAssemblies/{{sample}}/myloasm/assembly.fasta",
        contigs_info=f"{OUTDIR}/01-LongAssemblies/{{sample}}/myloasm/contigs_info.tsv",
        graph=f"{OUTDIR}/01-LongAssemblies/{{sample}}/myloasm/assembly_graph.gfa",
        alternate=f"{OUTDIR}/01-LongAssemblies/{{sample}}/myloasm/assembly_alternate.fa",
    params:
        tmpdir=f"{OUTDIR}/01-LongAssemblies/{{sample}}/myloasm/tmp",
        extra_flags=myloasm_extra_flags(),
    threads: config["threads"]
    conda:
        config["main_env"]
    log:
        f"{OUTDIR}/01-LongAssemblies/{{sample}}/myloasm/myloasm.log",
    shell:
        """
        mkdir -p {params.tmpdir}
        myloasm {input.reads} -o {params.tmpdir} -t {threads} --clean-dir {params.extra_flags} > {log} 2>&1
        cp {params.tmpdir}/assembly_primary.fa {output.assembly}
        {{ echo -e "seq_name\\tlength\\tcircular"; seqtk comp {output.assembly} \
            | awk -F"\\t" 'BEGIN{{OFS="\\t"}}{{circ=($1 ~ /circular-yes/) ? "Y" : "N"; print $1,$2,circ}}'; \
            }} > {output.contigs_info}
        mv {params.tmpdir}/final_contig_graph.gfa {output.graph}
        if [ -f {params.tmpdir}/alternate_assemblies/assembly_alternate.fa ]; then
            mv {params.tmpdir}/alternate_assemblies/assembly_alternate.fa {output.alternate}
        else
            touch {output.alternate}
        fi
        rm -rf {params.tmpdir}
        """
