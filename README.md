# plot-protein-variant

## Overview

Nextflow pipeline to fetch gene domain/transcript info and generate per-gene genomic/exonic/protein effect plots per gene.

## Requirements

- Nextflow version 22.10.0 or higher
- Singularity or Docker (or a compatible container runtime)

## Setup

- **Setup dependencies**:

```bash
make install
```

- **Prepare a gene-list**:

```text
GENE-1
GENE-2
```

- **Configure your reference and parameters** in `nextflow.config` or simply provide the parameters from the command line:

```groovy
params {
    // required
    gene_list         = ""  // one gene symbol per line
    outdir            = ""
    genomic           = ""  // coding genomic variants (delimiter: ;)
    gnomAD            = ""  // gnomAD variant frequencies (delimiter: ,)
    uk_Biobank        = ""  // control cohort variants (delimiter: \t)
    colour_phenotypes = ""  // phenotype -> colour mapping (XLSX)

    // optional - omit or leave blank to skip; the plot adapts to missing tracks
    genomic_noncoding = ""  // non-coding genomic variants (delimiter: ;)
    exonic            = ""  // exonic variants (delimiter: ;)
    sv                = ""  // structural variants (delimiter: ;)
    encode_file       = ""  // ENCODE functional elements (BED/GFF/GTF)
    refseq_file       = ""  // RefSeq annotation track (BED/GFF/GTF)
}
```

Any optional param left blank falls back to a placeholder input inside the pipeline, and the corresponding track/legend is simply left out of the plot for that gene.

Example usage:

```sh
nextflow run main.nf \
  --gene_list path/to/gene_list.txt \
  --genomic path/to/coding_variants.csv \
  --exonic path/to/exonic.csv \
  --gnomAD path/to/gnomAD.csv \
  --uk_Biobank path/to/Control_Cohort.tsv \
  --colour_phenotypes path/to/Color_Code.xlsx \
  --outdir results \
  -profile docker \
  -resume
```

## Output Directory Structure

```text
../results/
├── <GENE-1>
│   ├── <GENE-1>_domain.tsv
│   ├── <GENE-1>_transcript.tsv
│   ├── <GENE-1>_exonic_table.csv
│   ├── <GENE-1>_genomic_table.csv
│   ├── <GENE-1>_genomic_plot.png
│   ├── FETCH_GENE_INFO_log/.command.log
│   └── PLOT_VARIANTS_log/.command.log
├── <GENE-2>
│   ├── <GENE-2>_domain.tsv
│   ├── <GENE-2>_transcript.tsv
│   ├── <GENE-2>_exonic_table.csv
│   ├── <GENE-2>_genomic_table.csv
│   ├── <GENE-2>_genomic_plot.png
│   ├── FETCH_GENE_INFO_log/.command.log
│   └── PLOT_VARIANTS_log/.command.log
└── pipeline_info/
    ├── execution_timeline_<timestamp>.html
    ├── execution_report_<timestamp>.html
    ├── execution_trace_<timestamp>.txt
    └── pipeline_dag_<timestamp>.html
```

`pipeline_info/` is written once per run and covers the whole gene list (timeline, resource report, trace, DAG), not per-gene.

### Per-gene failure handling

Genes are processed independently, and one gene's failure never stops the run for the rest of the list:

- If a gene has no UniProt/Ensembl hit, `FETCH_GENE_INFO` is skipped for it and it produces no output.
- If a gene has zero variants passing the filtering thresholds, `<GENE>_genomic_table.csv` / `<GENE>_exonic_table.csv` are written as placeholder files containing `No variants found!` instead of a real table, and no plot is produced.

If a gene's output looks wrong or missing, check `<outdir>/<gene>/FETCH_GENE_INFO_log/.command.log` and `<outdir>/<gene>/PLOT_VARIANTS_log/.command.log` first.

## Input File Schema

A file's extension does not indicate how it's parsed — match the actual format and delimiter below, not just the column names. For example, `--exonic` is parsed as plaintext, not Excel, and `--UK_Biobank` is tab-delimited, not comma-delimited.

This covers only the params you pass to `nextflow run main.nf`. The pipeline also derives a couple of gene-specific inputs (protein domain coordinates, canonical transcript ID) on its own per gene, so there's nothing extra to supply for those.

| Param | Format | Delimiter | Required columns |
| --- | --- | --- | --- |
| `--genomic` | plaintext | `;` | `Sample-ID`, `Phenotype_complete`, `Chr:Pos`, `RefSeq Genes 110, NCBI`, `Effect (Combined)`, `HGVS c. (Clinically Relevant)`, `HGVS p. (Clinically Relevant)`, `gnomAD Genomes Variant Frequencies 4.0 v2, BROAD`, `PHRED-Score`, `REVEL-Score`, `Alt Allele Counts (AC)` |
| `--exonic` | plaintext | `;` | same columns as `--genomic` |
| `--genomic_noncoding` | plaintext | `;` | `Sample-ID`, `Phenotype_complete`, `Chr:Pos`, `RefSeq Genes 110, NCBI`, `Effect (Combined)`, `HGVS c. (Clinically Relevant)`, `Alt Allele Counts (AC)` (no gnomAD/PHRED/REVEL/AC filtering is applied to this file) |
| `--sv` | plaintext | `;` | `Sample-ID`, `Phenotype_complete`, `Chr:Pos`, `RefSeq Genes 110, NCBI`, `Effect (Combined)`, `HGVS g. (Clinically Relevant)`, `Alt Allele Counts (AC)` |
| `--gnomAD` | plaintext | `,` | `chrom`, `pos`, `AF_nfe` |
| `--UK_Biobank` | plaintext | `\t` | `SYMBOL`, `old_gnomAD_AF` |
| `--colour_phenotypes` | xlsx | | `Phenotype_complete`, `Colorcode` |
| `--encode_file` | genome annotation track (BED/GFF/GTF, format auto-detected from the file extension) | | standard track fields (chrom/start/end); a 4th metadata column or `name` column is used as the element label if present |
| `--refseq_file` | genome annotation track (BED/GFF/GTF, format auto-detected from the file extension) | | same as `--encode_file` |

`--genomic`, `--exonic`, `--genomic_noncoding`, and `--sv` also expect a comma as the decimal separator in numeric fields (e.g. `20,5`, not `20.5`).

`--genomic` and `--exonic` rows are additionally matched to the gene by a regex on `RefSeq Genes 110, NCBI`, and a variant is kept in the plot only if it passes all of:

- `gnomAD Genomes Variant Frequencies 4.0 v2, BROAD` &le; 0.001 (or missing)
- `PHRED-Score` &ge; 20 **or** `REVEL-Score` &ge; 0.5
- `Alt Allele Counts (AC)` &le; 4

These thresholds are fixed inside the pipeline, not exposed as CLI options. If a gene produces no plot, check its input rows against these cutoffs before assuming the pipeline is broken.

## Customization

- Adjust resource requirements and containers in `nextflow.config` and `conf/base.config`.
- Add or remove modules as needed for your workflow.

## Components

- Containers
  - `ghcr.io/maurya-anand/plot-protein-variant`
  - `python:3.11`

- Reference annotation
  - Ensembl `homo_sapiens`, genome build GRCh38, release 115

- R packages
  - BiocManager
  - colorspace
  - cowplot
  - dplyr
  - data.table
  - ensembldb
  - GenomicRanges
  - ggplot2
  - ggrepel
  - optparse
  - ragg
  - readr
  - readxl
  - rtracklayer
  - stringr
  - tidyr
