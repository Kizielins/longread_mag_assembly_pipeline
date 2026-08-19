# Pipeline

This strain-focused pipeline assembles MAGs from long-read ONT metagenomic data.

Two coupled outputs from one assembly: **Track A** is the species-level MAG catalog (assembly through binning/refinement/QC/taxonomy); **Track B** is strain resolution ([Strainy](https://github.com/katerinakazantseva/strainy), opt-in via `run_strain_phasing`, works with whichever assembler Track A uses).

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

## Installation instructions

```
mamba env create -f environment.yml            ## main env (flye, metabat2, medaka, DAS_Tool, SemiBin2, ...)
mamba env create -f environment.checkm2.yml    ## CheckM2, isolated due to a conflicting python pin
conda activate pipeline
```

## Database download

`download_databases.sh` fetches everything the pipeline needs in one pass: the
CheckM2 database, the GUNC database, the [GTDB-Tk](https://gtdb.ecogenomic.org/downloads)
reference package (pinned to **release226 / R226** for reproducibility), and
hostile's human host-reference index. Total size is 150+ GB, so make sure you
have the disk space before running it.

```
conda activate pipeline
bash download_databases.sh /path/to/checkm2/db /path/to/gunc/db /path/to/gtdbtk/db
```
Each path argument is optional (defaults live under `~/databases/`). The pipeline
reads all three paths straight from `config/config.yaml` (`checkm2_db`, `gunc_db`,
`gtdbtk_data_path`). Their defaults already match this script's own default locations, so a default-location download just
works; if you passed a custom path above, override the matching config key too
(directly in `config/config.yaml`, or `--config checkm2_db=... gunc_db=... gtdbtk_data_path=...`).

## Configuration

The `config/config.yaml` parameters are the pipeline's full set of run options:

## Running the pipeline

```
conda activate pipeline
snakemake -s workflow/Snakefile --configfile config/config.yaml \
    --config samples_sheet=samples.tsv outdir=pipeline-out --cores 16
```

`samples.tsv` is a TSV with a `sample` and `fastq` column, one row per sample:

```
sample	fastq
sample1	reads/sample1.fastq.gz
sample2	reads/sample2.fastq.gz
```

To compare assemblers in one run, pass a list instead of the single-item default, e.g. `--config 'assemblers=["metaflye","nanomdbg","myloasm"]'`.

Add `--use-conda` to let Snakemake manage each rule's conda env itself (needed for CheckM2 to resolve correctly).
