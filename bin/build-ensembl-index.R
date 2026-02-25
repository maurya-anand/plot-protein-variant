#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(optparse))


option_list <- list(
  make_option(c("-r", "--release"), type = "integer", default = 115,
              help = "Ensembl release [default: %default]"),
  make_option(c("-s", "--species"), type = "character", default = "homo_sapiens",
              help = "Species (FTP folder, snake_case) [default: %default]"),
  make_option(c("-g", "--genome"), type = "character", default = "GRCh38",
              help = "Genome build [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = getwd(),
              help = "Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

release <- opt$release
species <- opt$species
genome  <- opt$genome
outdir  <- opt$outdir

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ensembldb", quietly = TRUE)) {
  BiocManager::install("ensembldb", ask = FALSE, update = FALSE)
}

species_cap <- paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))

gtf_gz_name <- sprintf("%s.%s.%d.gtf.gz", species_cap, genome, release)
gtf_name    <- sub("\\.gz$", "", gtf_gz_name)
sqlite_name <- sprintf("%s.%s.%d.sqlite", species_cap, genome, release)

gtf_gz <- file.path(outdir, gtf_gz_name)
gtf    <- file.path(outdir, gtf_name)
sqlite <- file.path(outdir, sqlite_name)

url <- sprintf(
  "https://ftp.ensembl.org/pub/release-%d/gtf/%s/%s",
  release, species, gtf_gz_name
)

if (!file.exists(gtf)) {
  if (!file.exists(gtf_gz)) {
    message("Downloading GTF...")
    system(sprintf('curl -fL -o "%s" "%s"', gtf_gz, url), intern = FALSE)
  }
  message("Extracting GTF...")
  system(sprintf('gunzip -f "%s"', gtf_gz), intern = FALSE)
} else {
  message("GTF already exists. Skipping download/extraction.")
}

if (!file.exists(sqlite)) {
  message("Building EnsDb SQLite...")
  db_file_path <- ensembldb::ensDbFromGtf(
    gtf = gtf,
    outfile = sqlite,
    organism = species_cap,
    genomeVersion = genome,
    version = release
  )
  message(sprintf("Created: %s", db_file_path))
} else {
  message("SQLite already exists. Skipping database build.")
  db_file_path <- normalizePath(sqlite)
}

cat(db_file_path, "\n")