# Pipeline

This strain-focused pipeline assembles MAGs from long-read ONT metagenomic data.

Two coupled outputs from one assembly: **Track A** is the species-level MAG catalog (assembly through binning/refinement/QC/taxonomy); **Track B** is strain resolution ([Strainy](https://github.com/katerinakazantseva/strainy), opt-in via `run_strain_phasing`, works with whichever assembler Track A uses) - see [Strain resolution (Track B)](#strain-resolution-track-b) below.

Pipeline accepts Dorado's fastq output (`samtools fastq` first if you only have BAM), make sure `medaka_model` matches the Dorado basecalling model you used.

## Pipeline architecture

The pipeline itself is a [Snakemake](https://snakemake.readthedocs.io/) workflow, run with `snakemake` directly:

```
workflow/
├── Snakefile           entry point: includes the rule files below, defines the final target
├── rules/
│   ├── common.smk       shared helpers (dynamic bin-list lookup, config parsing)
│   ├── qc.smk            NanoPlot + Chopper (length/quality filter) + hostile (host-read removal)
│   ├── assembly.smk      metaFlye | nanoMDBG | myloasm - picked via the `assemblers` config list
│   ├── strainphasing.smk Track B: Strainy strain phasing on the raw assembly graph - switch on via `run_strain_phasing`
│   ├── polishing.smk     single whole-assembly medaka polish
│   ├── binning.smk       read mapping -> coverage (jgi_summarize_bam_contig_depths), metabat2 + SemiBin2
│   │                     (+ optional COMEBin)
│   ├── refinement.smk    DAS_Tool
│   ├── binqc.smk         CheckM2 -> GUNC chimera filter (two-stage) + barrnap/tRNAscan-SE (MIMAG HQ)
│   └── taxonomy.smk      GTDB-Tk + final genome-summary table
└── scripts/              quality filtering, summary-table build
config/config.yaml        all pipeline parameters
```
