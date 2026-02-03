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
    gene_list = ""
    outdir = ""
    genomic = ""
    genomic_noncoding = ""
    exonic = ""
    sv = ""
    encode_file = ""
    refseq_file = ""
    gnomAD = ""
    uk_Biobank = ""
    colour_phenotypes = ""
}
```

Example usage:

```sh
nextflow run main.nf \
  --gene_list path/to/gene_list.txt \
  --genomic path/to/coding_variants.csv \
  --exonic path/to/exonic.xlsx \
  --gnomAD path/to/gnomAD.csv \
  --uk_Biobank path/to/Control_Cohort.csv \
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
│   └── <GENE-1>_genomic_plot.png
└── <GENE-2>
    ├── <GENE-2>_domain.tsv
    ├── <GENE-2>_transcript.tsv
    └── <GENE-2>_genomic_plot.png
```

## Customization

- Adjust resource requirements and containers in `nextflow.config` and `conf/base.config`.
- Add or remove modules as needed for your workflow.

## Components

- BiocManager
- biomaRt
- colorspace
- cowplot
- dplyr
- GenomicRanges
- ggplot2
- ggrepel
- optparse
- ragg
- readr
- readxl
- rtracklayer
- stringr
