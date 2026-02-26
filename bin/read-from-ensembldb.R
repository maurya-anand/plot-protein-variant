#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages(library(optparse))


option_list <- list(
  make_option(c("-d", "--db"), type = "character", default = "/data/Homo_sapiens.GRCh38.115.sqlite",
              help = "Ensembl release index file [default: %default]"),
  make_option(c("-t", "--transcript"), type = "character", default = NULL,
              help = "Transcript ID [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$transcript)) {
  stop("Error: Transcript ID (--transcript) is required.")
}

if (!file.exists(opt$db)) {
  stop(sprintf("Error: Database file not found at: %s", opt$db))
}

db_path <- opt$db
transcript_id <- opt$transcript

message(sprintf("Processing transcript: %s", transcript_id))

suppressPackageStartupMessages({
  library(ensembldb)
  library(GenomicRanges)
  library(data.table)
  library(DBI)
})

pkgs <- c(
  "ensembldb", "GenomicRanges", "data.table", "DBI"
)

for (pkg in pkgs) {
  message(paste("Loading:", pkg))
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

get_transcript_info <- function(sqlite_path, transcript_id) {
  edb <- EnsDb(sqlite_path)
  on.exit(DBI::dbDisconnect(dbconn(edb)), add = TRUE)
  
  # EnsDb stores ENST IDs as "tx_name"
  tx_filter <- TxNameFilter(transcript_id)
  
  # exons in transcript order; exon_rank is the "rank" you want
  exl <- exonsBy(edb, by="tx", filter=tx_filter)
  if (length(exl) == 0) stop("Transcript not found: ", transcript_id)
  ex <- exl[[1]]
  
  exon_positions <- data.table(
    ensembl_transcript_id = transcript_id,
    exon_chrom_start = start(ex),
    exon_chrom_end   = end(ex),
    rank             = mcols(ex)$exon_rank
  )
  setorder(exon_positions, rank)
  
  # CDS ranges for transcript
  cdsl <- cdsBy(edb, by="tx", filter=tx_filter)
  if (length(cdsl) == 0) {
    # non-coding transcript; return empty cds_df with correct columns
    cds_df <- data.table(
      ensembl_transcript_id = character(),
      cds_start = integer(),
      cds_end = integer(),
      exon_chrom_start = integer(),
      exon_chrom_end = integer(),
      rank = integer(),
      chromosome_name = character(),
      strand = integer()
    )
    return(list(exon_positions = as.data.frame(exon_positions),
                cds_df = as.data.frame(cds_df)))
  }
  cds <- cdsl[[1]]
  
  # Map each CDS segment to its exon rank via overlap with exons
  hits <- findOverlaps(cds, ex, type="within")
  if (length(hits) == 0) stop("No CDS segments mapped to exons for ", transcript_id)
  
  cds_dt <- data.table(
    cds_genomic_start = start(cds)[queryHits(hits)],
    cds_genomic_end   = end(cds)[queryHits(hits)],
    rank              = mcols(ex)$exon_rank[subjectHits(hits)],
    chromosome_name   = as.character(seqnames(cds)[queryHits(hits)]),
    strand_char       = as.character(strand(cds)[queryHits(hits)])
  )
  
  # Ensure transcript order: exon rank already reflects transcript order (including minus strand)
  setorder(cds_dt, rank, cds_genomic_start, cds_genomic_end)
  
  cds_dt[, seg_len := cds_genomic_end - cds_genomic_start + 1L]
  cds_dt[, cds_end := cumsum(seg_len)]
  cds_dt[, cds_start := cds_end - seg_len + 1L]
  
  # Join in exon genomic span (BioMart columns)
  cds_dt <- merge(
    cds_dt,
    exon_positions,
    by = "rank",
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  strand_num <- if (cds_dt$strand_char[1] == "+") 1L else -1L
  
  cds_df <- cds_dt[, .(
    ensembl_transcript_id = transcript_id,
    cds_start,
    cds_end,
    exon_chrom_start,
    exon_chrom_end,
    rank,
    chromosome_name,
    strand = strand_num
  )]
  
  list(
    exon_positions = as.data.frame(exon_positions),
    cds_df = as.data.frame(cds_df)
  )
}

# Example (offline):
# transcript_id <- "ENST00000261917"

res <- get_transcript_info(db_path, transcript_id)

if (is.null(res$exon_positions) || nrow(res$exon_positions) == 0) {
  message("Exon positions are empty for transcript: ", transcript_id)
} else {
  message(sprintf("Finished extracting exon positions for: %s", transcript_id))
}

if (is.null(res$cds_df) || nrow(res$cds_df) == 0) {
  message("CDS data frame is empty for transcript: ", transcript_id)
} else {
  message(sprintf("Finished extracting cds positions for: %s", transcript_id))
}

write.csv(res$exon_positions, "exon_positions_offline.csv", row.names=FALSE,quote = FALSE)
write.csv(res$cds_df, "cds_df_offline.csv", row.names=FALSE,quote = FALSE)

saveRDS(res$exon_positions, "exon_positions.rds")
saveRDS(res$cds_df, "cds_df.rds")
