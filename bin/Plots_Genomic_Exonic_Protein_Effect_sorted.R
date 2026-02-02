#!/usr/bin/env Rscript

pkgs <- c(
  "dplyr", "readxl", "rtracklayer", "GenomicRanges", "stringr",
  "ggplot2", "cowplot", "ggrepel", "biomaRt", "ragg", "readr", "optparse"
)
for (pkg in pkgs) {
  print(paste("Loading:", pkg))
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

###Dezimalschreibweise
options(scipen = 999)
options(timeout = 300)
##################################################################################

option_list = list(
    make_option(c("--gene_name"), type="character", default=NULL, help="", metavar="character"),
    make_option(c("--gene_domain"), type="character", default=NULL, help="", metavar="character"),
    make_option(c("--transcript_id"), type="character", default=NULL, help="", metavar="character"),
    make_option(c("--genomic"), type="character", default=NULL, help="", metavar="character"),
    make_option(c("--genomic_noncoding"), type="character", default=NULL, help="", metavar="character"), 
    make_option(c("--exonic"), type="character", default=NULL, help="", metavar="character"),
    make_option(c("--sv"), type="character", default=NULL, help="[default= %default]", metavar="character"),
    make_option(c("--encode_file"), type="character", default=NULL, help="[default= %default]", metavar="character"), 
    make_option(c("--refseq_file"), type="character", default=NULL, help="[default= %default]", metavar="character"),
    make_option(c("--gnomAD"), type="character", default=NULL, help="[default= %default]", metavar="character"),
    make_option(c("--UK_Biobank"), type="character", default=NULL, help="[default= %default]", metavar="character"),
    make_option(c("--colour_phenotypes"), type="character", default=NULL, help="[default= %default]", metavar="character") 
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# ###Individuell anpassbare Parameter
# sanitize_file_arg <- function(x) {
#   if (is.null(x) || !nzchar(x) || grepl("^NO_FILE_", x)) return(NULL)
#   return(x)
# }

# genomic           <- sanitize_file_arg(opt$genomic)
# genomic_noncoding <- sanitize_file_arg(opt$genomic_noncoding)
# exonic            <- sanitize_file_arg(opt$exonic)
# sv                <- sanitize_file_arg(opt$sv)
# encode_file       <- sanitize_file_arg(opt$encode_file)
# refseq_file       <- sanitize_file_arg(opt$refseq_file)
# gnomAD            <- sanitize_file_arg(opt$gnomAD)
# UK_Biobank        <- sanitize_file_arg(opt$UK_Biobank)
# colour_phenotypes <- sanitize_file_arg(opt$colour_phenotypes)
##Dateipfade einlesen
genomic           <- opt$genomic
genomic_noncoding <- opt$genomic_noncoding
exonic            <- opt$exonic
sv <- opt$sv
encode_file <- opt$encode_file
refseq_file <- opt$refseq_file

# # if (!is.null(opt$genomic_noncoding)){
# #   genomic_noncoding <- opt$genomic_noncoding
# # }
# # if (!is.null(opt$sv)){
# #   sv <- opt$sv
# # }
# # if (!is.null(opt$encode_file)){
# #   encode_file <- opt$encode_file
# # }
# # if (!is.null(opt$refseq_file)){
# #   refseq_file <- opt$refseq_file
# # }
gnomAD            <- opt$gnomAD
UK_Biobank        <- opt$UK_Biobank
colour_phenotypes <- opt$colour_phenotypes

##Filterkriterien_coding
gene_name <- opt$gene_name
af_max    <- 0.0001
cadd_min  <- 25
revel_min <- 0.5
ac_max    <- 4

##Plot-Parameter
EXON_SCALE_FACTOR  <-  0.3   #Exonbreite - relativ zur echten genomischen Länge 
INTRON_GAP         <-  250   #Intronbreite
EXON_HALF_HEIGHT   <-  0.03  #Exon-Höhe - Gesamthöhe 2*Exon_Half_Height
STACK_DY           <-  0.05  #vertikaler Abstand zwischen gestapelten Varianten
X_STACK_TOL        <-  400   #Breite eines horizontalen "x-Clusters" in der nahe beieinanderliegende Varianten gestapelt werden - je größer, desto mehr Varianten werden gestapelt
Y_GENOMIC_NC       <-  0.60  #Höhe NC-Varianten Genom
Y_GENOMIC_LOF      <-  0.45  #Höhe LoF-Varianten Genom
Y_GENOMIC_MISS     <-  0.30  #Höhe Missense-Varianten Genom
Y_GENOMIC_OTHER    <-  0.15  #Höhe Other-Varianten Genom
Y_EXOM_OTHER       <- -0.15  #Höhe Other-Varianten Exom
Y_EXOM_MISS        <- -0.30  #Höhe Missense-Varianten Exom
Y_EXOM_LOF         <- -0.45  #Höhe LoF-Varianten Exom
LABEL_NUDGE_X      <-  0     #Horizontale Position der Labels vom Punkt
LABEL_ANGLE        <-  2     #Winkel der Label
TICK_LEN           <-  0.1   #Länge der Ticks coding nach oben/unten
SV_HALF_HEIGHT     <-  0.002 #Höhe der SV-Darstellung

Y_FUNC_RECT_BOTTOM <- -0.03  #Anpassung der Boxen für funktionelle Elemente/Nähe zum Gen
Y_FUNC_RECT_TOP    <-  0.03  #Anpassung der Boxen für funktionelle Elemente/Höhe der Boxen

GENE_LINE_WIDTH    <-  1.5   #Linien im Genplot
EXON_BORDER_WIDTH  <-  1.5   #Linien im Genplot
EXON_LABEL_SIZE    <-  5.0   #Schriftgröße Exon-Nr
DOMAIN_HALF_HEIGHT <- EXON_HALF_HEIGHT * 0.89 #Höhen der DOmänen (90% vom Exom)

Y_LABEL_POS  <- c(           #Position/Reihenfolge der Labels auf Y-Achse
  Y_GENOMIC_NC,
  Y_GENOMIC_LOF,
  Y_GENOMIC_MISS,
  Y_GENOMIC_OTHER,
  Y_EXOM_OTHER,
  Y_EXOM_MISS,
  Y_EXOM_LOF
)

Y_LABEL_TEXT <- c(           #Name der Labels auf Y-Achse
  "non-coding (GS)",
  "LoF (GS)",
  "Missense (GS)",
  "Other (GS)",
  "Other (ES)",
  "Missense (ES)",
  "LoF (ES)"
)

##Density-Plot
DENSITY_HEIGHT    <- 0.10    # Höhe des Density-Tracks
DENSITY_N         <- 1024    # Auflösung der Kurve
DENSITY_ADJUST    <- 1       # Glättung (größer = glatter)

##Ensembl-Daten
transcript_id <- opt$transcript_id #Individuell an das Gene of interest anpassen

domain_tbl <- read_tsv(opt$gene_domain, show_col_types = FALSE)
domain_input <- domain_tbl %>%
  dplyr::select(domain = name, aa_start = start, aa_end = end)
domain_names <- unique(domain_input$domain)
palette <- colorspace::sequential_hcl(length(domain_names), palette = "TealGrn")
domain_colors <- setNames(palette, domain_names)
domain_input$color <- domain_colors[domain_input$domain]

##Plottitel
plot_title <- paste0(
  gene_name, " - Variant display\n",
  "AF (gnomAD v3 & v4) < ", af_max,
  ", REVEL > ", revel_min,
  ", CADD > ", cadd_min
)

##nach dem Part Datein einlesen, vorberieten und filtern muss nichts mehr individuell angepast werden, außer man möchte die Schriftgröße oder so nochmal ändern
##################################################################################
##Helper fürs einlesen der Datein
read_table_auto <- function(path, na = c("", "NA"), ...) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("xlsx", "xls")) {
    return(readxl::read_excel(path = path, na = na, ...))
  }
  
  if (ext %in% c("csv", "tsv", "txt")) {
    delim <- if (ext == "tsv") "\t" else ","
    return(readr::read_delim(
      file = path,
      delim = delim,
      na = na,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "minimal",   #WICHTIG: Namen NICHT reparieren
      ...
    ))
  }
  
  stop("Unbekanntes Dateiformat: ", ext, " (", path, ")")
}

##################################################################################
###Daten einlesen, vorbereiten und filtern
###Wenn die Datei nicht existiert, wird sie auf Null gesetzt und das Skript geht trotzdem
###Wichtig: Pfad darf nicht auskommentiert sein!

##Genom-Daten vorbereiten und verkleinern - Datei einlesen, wichtige Spalten nummerisch machen, Datensatz verkleinern
##Spaltennamen entsprechend an Datei anpassen
if (file.exists(genomic)) {
  message("Genomic-Datei gefunden")
  genomic <- read_table_auto (genomic, na = c("", "NA"))
  genomic$`gnomAD Genomes Variant Frequencies 4.0 v2, BROAD` <- as.numeric(genomic$`gnomAD Genomes Variant Frequencies 4.0 v2, BROAD`)
  genomic$`PHRED-Score` <- as.numeric(genomic$`PHRED-Score`)
  genomic$`REVEL-Score` <- as.numeric(genomic$`REVEL-Score`)
  
  genomic_small <- genomic %>%
    dplyr::select(`Sample-ID`, 
                  `Phenotype_complete`, 
                  `Chr:Pos`, 
                  `RefSeq Genes 110, NCBI`,
                  `Effect (Combined)`,
                  `HGVS c. (Clinically Relevant)`,
                  `HGVS p. (Clinically Relevant)`,
                  `gnomAD Genomes Variant Frequencies 4.0 v2, BROAD`,
                  `PHRED-Score`,
                  `REVEL-Score`,
                  `Alt Allele Counts (AC)`)
} else {
  message("Genomic-Datei nicht gefunden")
  genomic_small <- NULL
}

##Nicht-kodierende Genom-Daten vorbereiten und verkleinern - Datei einlesen, wichtige Spalten nummerisch machen, Datensatz verkleinern
##SPaltennamen entsprechend an Datei anpassen
if (file.exists(genomic_noncoding)) {
  message("Noncoding-Datei gefunden")
  genomic_noncoding <- read_table_auto(genomic_noncoding, na = c("", "NA"))
  
  genomic_noncoding_small <- genomic_noncoding %>%
    dplyr::select(
      `Sample-ID`,
      `Phenotype_complete`,
      `Chr:Pos`,
      `RefSeq Genes 110, NCBI`,
      `Effect (Combined)`,
      `HGVS c. (Clinically Relevant)`,
      `Alt Allele Counts (AC)`
    ) %>% 
    filter(grepl(paste0("\\b", gene_name, "\\b"), `RefSeq Genes 110, NCBI`))
} else {
  message("Noncoding-Datei nicht gefunden")
  genomic_noncoding_small <- NULL
}


##SV-Daten vorbereiten und verkleinern - Datei einlesen, wichtige Spalten nummerisch machen, Datensatz verkleinern
##Spaltennamen entsprechend an Datei anpassen
if (file.exists(sv)) {
  message("SV-Datei gefunden")
  sv_raw <- read_table_auto(sv, na = c("", "NA"))
  
  sv_small <- sv_raw %>%
    dplyr::select(
      `Sample-ID`,
      `Phenotype_complete`,
      `Chr:Pos`,
      `RefSeq Genes 110, NCBI`,
      `Effect (Combined)`,
      `HGVS g. (Clinically Relevant)`,
      `Alt Allele Counts (AC)`
    ) %>% 
    filter(grepl(paste0("\\b", gene_name, "\\b"), `RefSeq Genes 110, NCBI`))
  
} else {
  message("SV-Datei nicht gefunden")
  sv_small <- NULL
}

##Exom-Daten vorbereiten und verkleinern - Datei einlesen, wichtige Spalten nummerisch machen, Datensatz verkleinern
##Spaltennamen entsprechend an Datei anpassen
if (file.exists(exonic)) {
  message("Exonic-Datei gefunden")
  exonic <- read_table_auto(exonic, na = c("", "NA"))
  exonic$`gnomAD Genomes Variant Frequencies 4.0 v2, BROAD` <- as.numeric(exonic$`gnomAD Genomes Variant Frequencies 4.0 v2, BROAD`)
  exonic$`PHRED-Score` <- as.numeric(exonic$`PHRED-Score`)
  exonic$`REVEL-Score` <- as.numeric(exonic$`REVEL-Score`)
  
  exonic_small <- exonic %>%
    dplyr::select(`Sample-ID`, 
                  `Phenotype_complete`, 
                  `Chr:Pos`, 
                  `RefSeq Genes 110, NCBI`,
                  `Effect (Combined)`,
                  `HGVS c. (Clinically Relevant)`,
                  `HGVS p. (Clinically Relevant)`,
                  `gnomAD Genomes Variant Frequencies 4.0 v2, BROAD`,
                  `PHRED-Score`,
                  `REVEL-Score`,
                  `Alt Allele Counts (AC)`)
} else {
  message("Exonic-Datei nicht gefunden")
  exonic_small <- NULL
}

##Genom-Daten coding filtern
if (!is.null(genomic_small) && nrow(genomic_small) > 0) {
  message("Genomic_small-Datei gefunden")
  genomic_coding_filtered <- genomic_small %>%
    filter(
      grepl(paste0("\\b", gene_name, "\\b"), `RefSeq Genes 110, NCBI`),               
      (is.na(`gnomAD Genomes Variant Frequencies 4.0 v2, BROAD`) |
         `gnomAD Genomes Variant Frequencies 4.0 v2, BROAD` <= af_max),  
      #(`PHRED-Score` >= cadd_min | `REVEL-Score` >= revel_min),      #Filterung REVEL oder CADD
      `PHRED-Score` >= cadd_min,
      `REVEL-Score` >= revel_min,
      `Alt Allele Counts (AC)` <= ac_max)
} else {
  genomic_coding_filtered <- NULL 
}

##Exom-Daten coding filtern 
if (!is.null(exonic_small) && nrow(exonic_small) > 0) {
  message("Exonic_small-Datei gefunden")
  exonic_coding_filtered <- exonic_small %>%   
    filter(
      grepl(paste0("\\b", gene_name, "\\b"), `RefSeq Genes 110, NCBI`),
      (is.na(`gnomAD Genomes Variant Frequencies 4.0 v2, BROAD`) |
         `gnomAD Genomes Variant Frequencies 4.0 v2, BROAD` <= af_max),
      #(`PHRED-Score` >= cadd_min | `REVEL-Score` >= revel_min),      #Filterung REVEL oder CADD
      `PHRED-Score` >= cadd_min,
      `REVEL-Score` >= revel_min,
      `Alt Allele Counts (AC)` <= ac_max)
} else {
  exonic_coding_filtered <- NULL
}

##UK Biobank einlesen + filtern
UKB_small <- NULL

if (file.exists(UK_Biobank)) {
  message("UK Biobank-Datei gefunden")
  
  UKB_raw <- readr::read_tsv(UK_Biobank, show_col_types = FALSE, progress = FALSE)
  
  UKB_small <- UKB_raw %>%
    mutate(
      SYMBOL = toupper(stringr::str_trim(as.character(SYMBOL))),
      old_gnomAD_AF = as.numeric(gsub(",", ".", as.character(old_gnomAD_AF)))
      #n_het_controls = suppressWarnings(as.numeric(n_het_controls))
    ) %>%
    filter(
      SYMBOL == toupper(gene_name),
      is.na(old_gnomAD_AF) | old_gnomAD_AF < 0.001
      #is.na(n_het_controls) | n_het_controls > 1
    )
  
  nrow(UKB_small)
  
} else {
  message("UK Biobank-Datei nicht gefunden")
}

##GnomAD einlesen + filtern
gnomAD_small <- NULL

if (!is.null(gnomAD) && file.exists(gnomAD)) {
  message("gnomAD-Datei gefunden")
  
  gnomAD_raw <- readr::read_csv(gnomAD, show_col_types = FALSE, progress = FALSE)
  
  gnomAD_small <- gnomAD_raw %>%
    mutate(
      chrom   = gsub("^chr", "", as.character(.data$chrom), ignore.case = TRUE),
      pos     = as.numeric(.data$pos),
      var_key = paste0(.data$chrom, ":", .data$pos),
      AF_nfe  = as.numeric(gsub(",", ".", as.character(.data$AF_nfe)))
      ) %>%
   filter(
    !is.na(.data$pos),
    is.na(.data$AF_nfe) | .data$AF_nfe < 0.001
     )
  
  message("gnomAD_small n = ", nrow(gnomAD_small))
  } else {
  message("gnomAD Datei nicht gefunden")
  }

##non-coding und SV-Datei wird gefiltert eingelesen
##es gibt einige Händische Filterschritte die über verschiedene Portale durchgeführt werden müssen

##################################################################################
###Plotgrundlagen-Vorbereitung

##Ensembl / Exon-Informationen laden
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", verbose = TRUE)

exon_positions <- getBM(
  attributes = c("ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end", "rank"),
  filters    = "ensembl_transcript_id",
  values     = transcript_id,
  mart       = ensembl
)

if (nrow(exon_positions) == 0) {
  stop("Für das Transkript ", transcript_id, " wurden keine Exon-Informationen gefunden.")
}

##Exons skalieren
##Exons haben skalierte Längen und zwischen denen liegen die Introns
exon_df <- exon_positions %>%
  dplyr::rename(exon_start = exon_chrom_start,
                exon_end   = exon_chrom_end,
                exon_id    = rank) %>%   
  arrange(exon_start) %>%                
  mutate(
    exon_width   = exon_end - exon_start + 1,
    scaled_width = pmax(1, round(exon_width * EXON_SCALE_FACTOR))
  )

n_exons <- nrow(exon_df)

scaled_start <- numeric(n_exons)
scaled_end   <- numeric(n_exons)

cur <- 1
for (i in seq_len(n_exons)) {
  scaled_start[i] <- cur
  scaled_end[i]   <- cur + exon_df$scaled_width[i] - 1
  if (i < n_exons) cur <- scaled_end[i] + 1 + INTRON_GAP
}

exon_df$scaled_start <- scaled_start
exon_df$scaled_end   <- scaled_end

exon_plot_df <- exon_df

TX_MIN <- min(exon_plot_df$scaled_start, na.rm = TRUE)
TX_MAX <- max(exon_plot_df$scaled_end,   na.rm = TRUE)

if (!is.finite(TX_MIN) || !is.finite(TX_MAX) || TX_MIN >= TX_MAX) {
  stop("Ungültige Exon-Koordinaten: TX_MIN/TX_MAX sind nicht sinnvoll.")
}

##Mapping: genomische Position
##genomische Basenpositionen werden auf Exon-Intron-Plot-Achse projiziert
map_genomic_to_scaled <- function(pos_vec) {
  sapply(pos_vec, function(pos) {
    hit <- which(pos >= exon_df$exon_start & pos <= exon_df$exon_end)
    if (length(hit) == 1) {
      j <- hit
      frac <- (pos - exon_df$exon_start[j]) /
        (exon_df$exon_end[j] - exon_df$exon_start[j] + 1)
      return(exon_df$scaled_start[j] +
               frac * (exon_df$scaled_end[j] - exon_df$scaled_start[j]))
    }
    left  <- max(which(pos > exon_df$exon_end), 0)
    right <- min(which(pos < exon_df$exon_start), n_exons + 1)
    
    if (left > 0 && right <= n_exons) {
      intron_start <- exon_df$scaled_end[left] + 1
      intron_end   <- exon_df$scaled_start[right] - 1
      frac <- (pos - exon_df$exon_end[left]) /
        (exon_df$exon_start[right] - exon_df$exon_end[left])
      return(intron_start + frac * (intron_end - intron_start))
    }
    
    if (pos < min(exon_df$exon_start))
      return(exon_df$scaled_start[1] - INTRON_GAP/2)
    
    return(exon_df$scaled_end[n_exons] + INTRON_GAP/2)
  })
}

##UKB Positionen -> scaled_x für Density
UKB_pos <- UKB_small %>%
  mutate(
    Location = as.character(Location),
    pos = as.numeric(stringr::str_extract(Location, "(?<=:)\\d+"))
  ) %>%
  filter(!is.na(pos)) %>%
  pull(pos)

UKB_scaled_x <- as.numeric(
  unlist(map_genomic_to_scaled(UKB_pos), use.names = FALSE)
)

UKB_scaled_x <- UKB_scaled_x[is.finite(UKB_scaled_x)]

length(UKB_scaled_x)
summary(UKB_scaled_x)

##gnomAD Positionen -> scaled_x für Density
gnomAD_pos <- gnomAD_small %>%
  filter(!is.na(pos)) %>%
  pull(pos)

gnomAD_scaled_x <- as.numeric(
  unlist(map_genomic_to_scaled(gnomAD_pos), use.names = FALSE)
)

gnomAD_scaled_x <- gnomAD_scaled_x[is.finite(gnomAD_scaled_x)]

length(gnomAD_scaled_x)
summary(gnomAD_scaled_x)

##Domänen (AA-Positionen → genomisch → Plot-Koordinaten)
##CDS-Informationen für das aktuelle Transkript holen
cds_df <- getBM(
  attributes = c("ensembl_transcript_id",
                 "cds_start", "cds_end", 
                 "exon_chrom_start", "exon_chrom_end",
                 "rank", "chromosome_name", "strand"),
  filters = "ensembl_transcript_id",
  values  = transcript_id,
  mart    = ensembl
)

cds_df <- cds_df %>% arrange(cds_start)
strand <- unique(cds_df$strand)
if (length(strand) != 1) stop("Fehler: mehrere Strands gefunden.")
strand <- strand[1]

message("Strand: ", strand)

##Funktion: AA-Range wird in CDS-Koordinaten übersetzt und dann in genomische Koordinaten 
##Idee: korrekte Position der Domänen wird später im Exon angezeigt
aa_to_genomic <- function(aa_start, aa_end, cds_df, strand) {
  
  genomic_coords <- list()
  cds_start_pos <- (aa_start - 1) * 3 + 1
  cds_end_pos   <- aa_end * 3
  cds_cursor <- 0
  
  for (i in seq_len(nrow(cds_df))) {
    
    exon_len <- cds_df$cds_end[i] - cds_df$cds_start[i] + 1
    
    #Exon ohne CDS: überspringen für Domänen
    if (is.na(cds_df$cds_start[i]) || is.na(cds_df$cds_end[i])) {
      next
    }
    
    cds_s <- cds_cursor + 1
    cds_e <- cds_cursor + exon_len
    
    overlap_start <- max(cds_start_pos, cds_s)
    overlap_end   <- min(cds_end_pos, cds_e)
    
    if (!is.na(overlap_start) && overlap_start <= overlap_end) {
      
      offset <- overlap_start - cds_s
      
      if (strand == "+") {
        genomic_start <- cds_df$exon_chrom_start[i] + offset
        genomic_end   <- genomic_start + (overlap_end - overlap_start)
      } else {
        genomic_end   <- cds_df$exon_chrom_end[i] - offset
        genomic_start <- genomic_end - (overlap_end - overlap_start)
      }
      
      genomic_coords[[length(genomic_coords) + 1]] <-
        c(genomic_start, genomic_end)
    }
    
    cds_cursor <- cds_cursor + exon_len
  }
  
  genomic_coords
}

##Für jede Domäne genomische + skalierte Koordinaten berechnen
domain_list <- list()

for (i in seq_len(nrow(domain_input))) {
  d <- domain_input[i, ]
  
  genomic_ranges <- aa_to_genomic(
    aa_start = d$aa_start,
    aa_end   = d$aa_end,
    cds_df   = cds_df,
    strand   = strand
  )
  
  if (is.null(genomic_ranges)) next
  
  genomic_ranges_df <- do.call(rbind, lapply(genomic_ranges, function(x) {
    data.frame(genomic_start = x[1], genomic_end = x[2])
  }))
  
  genomic_ranges_df$scaled_start <- map_genomic_to_scaled(genomic_ranges_df$genomic_start)
  genomic_ranges_df$scaled_end   <- map_genomic_to_scaled(genomic_ranges_df$genomic_end)
  
  genomic_ranges_df$domain <- d$domain
  genomic_ranges_df$color  <- d$color
  
  domain_list[[i]] <- genomic_ranges_df
}

##alles zu einem Dataframe zusammenfassen
##domain_df enthält genomische & skalierte Grenzen, Domänenname, Farben
if (length(domain_list) > 0) {
  domain_df <- dplyr::bind_rows(domain_list)
  domain_df$scaled_start <- as.numeric(domain_df$scaled_start)
  domain_df$scaled_end   <- as.numeric(domain_df$scaled_end)
} else {
  domain_df <- NULL
}

##Phenotyp-Farben laden aus einer Excel-Datei
##Datei enthält Phänotypen und den Hexadezimalcode für die Farbe
if (file.exists(colour_phenotypes)) {
  colmap_df <- read_excel(colour_phenotypes)
  colname_color <- intersect(("Colorcode"), names(colmap_df))[1]
  colname_pheno <- intersect(("Phenotype_complete"), names(colmap_df))[1]
  phenotype_colors <- setNames(colmap_df[[colname_color]],
                               colmap_df[[colname_pheno]])
} else {
  phenotype_colors <- NULL
}

##Hilfsfunktion für saubere Strings
##entfernt Steuerzeichen, konvertiert nach UTF-8, trimmt Whitespace
clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("[[:cntrl:]]", " ", x)
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- gsub(" +", " ", x)
  trimws(x)
}

##################################################################################
###Variantendaten-Vorbereitung

##Varianten-Dataframe vorbereiten
prepare_variant_df <- function(df) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df <- df %>%
    mutate(`Chr:Pos`= trimws(as.character(`Chr:Pos`))) %>%
    tidyr::separate(
      `Chr:Pos`, 
      into = c("chrom_raw","pos"),
      sep = ":", 
      convert = TRUE
    ) %>%
    mutate(
      chrom = gsub("^chr", "", chrom_raw, ignore.case = TRUE),
      var_key = paste0(chrom, ":", pos),
      scaled_x = map_genomic_to_scaled(pos))
  
  if ("Phenotype_complete" %in% names(df)) {
    df$Phenotype_complete <- clean_text(df$Phenotype_complete)
  }
  
   col_hgvs_p <- intersect("HGVS p. (Clinically Relevant)", names(df))[1]
   col_hgvs_c <- intersect("HGVS c. (Clinically Relevant)", names(df))[1]
  
   if (!is.na(col_hgvs_p)) {
     df$label_raw <- clean_text(df[[col_hgvs_p]])
     df$label     <- stringr::str_extract(df$label_raw, "p\\.[^,;\\s]+")
   } else if (!is.na(col_hgvs_c)) {
     df$label_raw <- clean_text(df[[col_hgvs_c]])
     df$label     <- stringr::str_extract(df$label_raw, "c\\.[^,;\\s]+")
   } else {
     df$label <- NA_character_
   }
  
  df$label <- clean_text(df$label)
  df$label[df$label == ""] <- NA_character_
  
  if ("Effect (Combined)" %in% names(df)) {
    eff     <- as.character(df$`Effect (Combined)`)
    eff_low <- tolower(eff)
    
    df$Effect_simple <- dplyr::case_when(
      eff_low == "missense" ~ "Missense",
      eff_low == "lof"      ~ "LOF",
      TRUE                  ~ "Other/Unknown"
    )
  } else {
    df$Effect_simple <- "Other/Unknown"
  }
  
  df$Effect_simple <- factor(df$Effect_simple,
                             levels = c("LOF", "Missense", "Other/Unknown"))
  
  df$is_sv        <- FALSE
  df$scaled_start <- NA_real_
  df$scaled_end   <- NA_real_
  
  return(df)
}

##SV-Daten vorbereiten (Block statt Punkt)
prepare_sv_df <- function(df) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df <- df %>%
    mutate(`Chr:Pos` = trimws(as.character(`Chr:Pos`))) %>%
    mutate(
      chrom_raw = as.integer(gsub("^chr", "", sub(":.*", "", `Chr:Pos`))),
      range_str = sub(".*:", "", `Chr:Pos`),
      start     = as.numeric(sub("-.*", "", range_str)),
      end       = as.numeric(sub(".*-", "", range_str)),
      pos_mid   = floor((start + end) / 2),
      chrom     = as.character(chrom_raw),
      var_key   = paste0(chrom, ":", pos_mid),
      scaled_x      = map_genomic_to_scaled(pos_mid),
      scaled_start  = map_genomic_to_scaled(start),
      scaled_end    = map_genomic_to_scaled(end)
    )
  
  if ("Phenotype_complete" %in% names(df)) {
    df$Phenotype_complete <- clean_text(df$Phenotype_complete)
  }
  
  col_hgvs_g <- intersect("HGVS g. (Clinically Relevant)", names(df))[1]
  
  if (!is.null(col_hgvs_g) && !is.na(col_hgvs_g)) {
    df$hgvs_g_raw   <- clean_text(df[[col_hgvs_g]])
    df$hgvs_g_label <- stringr::str_extract(df$hgvs_g_raw, "g\\.[^,;\\s]+")
  } else {
    df$hgvs_g_label <- NA_character_
  }
  
  df$label <- clean_text(df$hgvs_g_label)
  df$label[df$label == ""] <- NA_character_
  
  df$Effect_simple <- factor(
    "Other/Unknown",
    levels = c("Non-coding", "LOF", "Missense", "Other/Unknown")
  )
  
  df$is_sv <- TRUE
  
  df
}

##"Master"-Varianten-Tabellen vorbereiten
##es wird immer geprüft, ob die Datein erstellt wurden, wenn sie nicht existieren, werden sie auf Null gesetzt
vars_genomic_coding <- if (!is.null(genomic_coding_filtered) && nrow(genomic_coding_filtered) > 0) {
  prepare_variant_df(genomic_coding_filtered)
} else {
  NULL
}

vars_genomic_noncoding <- if (!is.null(genomic_noncoding_small) && nrow(genomic_noncoding_small) > 0) {
  tmp <- prepare_variant_df(genomic_noncoding_small)
  tmp$Effect_simple <- "Non-coding"
  tmp$Effect_simple <- factor(
    tmp$Effect_simple,
    levels = c("LOF", "Missense", "Non-coding", "Other/Unknown")
  )
  
  tmp
} else {
  NULL
}

vars_exonic <- if (!is.null(exonic_coding_filtered) && nrow(exonic_coding_filtered) > 0) {
  prepare_variant_df(exonic_coding_filtered)
} else {
  NULL
}

vars_sv <- if (!is.null(sv_small) && nrow(sv_small) > 0) {
  prepare_sv_df(sv_small)
} else {
  NULL
}

##Alles, was im Genom und Exom vorkommt, bleibt NUR im Genom
if(!is.null(vars_genomic_coding) && !is.null(vars_exonic)) {
  common_keys <- intersect(vars_genomic_coding$var_key,
                           vars_exonic$var_key)
  vars_exonic_unique <- vars_exonic %>%
    dplyr::filter(!var_key %in% common_keys)
} else {
  common_keys <- character(0)
  vars_exonic_unique <- vars_exonic
}

##coding und non-coding Varianten (SNVs und SVs) zu einer Genom-Liste zusammenbauen
variant_list <- list(
  vars_genomic_coding,
  vars_genomic_noncoding,
  vars_sv
)

variant_list <- Filter(Negate(is.null), variant_list)

variant_list <- lapply(variant_list, function(df) {
  if ("Alt Allele Counts (AC)" %in% names(df)) {
    df[["Alt Allele Counts (AC)"]] <- as.numeric(df[["Alt Allele Counts (AC)"]])
  }
  df
})

if (length(variant_list) == 0) {
  vars_genomic_all <- NULL
} else {
  vars_genomic_all <- dplyr::bind_rows(variant_list) %>%
    mutate(
      Effect_simple = as.character(Effect_simple),
      Effect_simple = dplyr::case_when(
        Effect_simple %in% c("LOF", "Missense", "Non-coding", "Other/Unknown") ~ Effect_simple,
        TRUE ~ "Other/Unknown"
      ),
      Effect_simple = factor(
        Effect_simple,
        levels = c("Non-coding", "LOF", "Missense", "Other/Unknown")
      )
    )
}

##################################################################################
###Datenaufbereitung für den Plot

##Genom: Coding und Non-coding trennen
vars_genomic_plot       <- NULL
vars_genomic_nc_plot    <- NULL
vars_exomic_plot        <- NULL
vars_genomic_nc_raw     <- NULL
vars_genomic_coding_raw <- NULL

if (!is.null(vars_genomic_all) && nrow(vars_genomic_all) > 0) {
  vars_genomic_nc_raw <- vars_genomic_all %>%
    dplyr::filter(Effect_simple == "Non-coding")
  
  vars_genomic_coding_raw <- vars_genomic_all %>%
    dplyr::filter(Effect_simple != "Non-coding")
}

vars_exonic_all <- vars_exonic_unique

##GENOM: pro Variante zusammenfassen
vars_genomic_plot <- NULL

if (!is.null(vars_genomic_all) && nrow(vars_genomic_all) > 0) {
  vars_genomic_plot <- vars_genomic_all %>%
    group_by(var_key) %>%
    summarise(
      scaled_x           = mean(scaled_x, na.rm = TRUE),
      scaled_start       = dplyr::first(scaled_start),   
      scaled_end         = dplyr::first(scaled_end),     
      Phenotype_complete = dplyr::first(Phenotype_complete),
      Effect_simple      = dplyr::first(Effect_simple),
      label              = dplyr::first(label),
      is_sv              = any(is_sv, na.rm = TRUE),      
      .groups = "drop"
    ) %>%
    mutate(
      in_both = var_key %in% common_keys
    )
}

##Non-coding Varianten für funktionellen Track
if (!is.null(vars_genomic_nc_raw) && nrow(vars_genomic_nc_raw) > 0) {
  vars_genomic_nc_plot <- vars_genomic_nc_raw %>%
    group_by(var_key) %>%
    summarise(
      scaled_x           = mean(scaled_x, na.rm = TRUE),
      Phenotype_complete = dplyr::first(Phenotype_complete),
      Effect_simple      = dplyr::first(Effect_simple),
      label              = dplyr::first(label),
      .groups = "drop"
    )
}

##EXOM-Varianten aggregieren
if (!is.null(vars_exonic_all) && nrow(vars_exonic_all) > 0) {
  vars_exomic_plot <- vars_exonic_all %>%
    group_by(var_key) %>%
    summarise(
      scaled_x           = mean(scaled_x, na.rm = TRUE),
      Phenotype_complete = dplyr::first(Phenotype_complete),
      Effect_simple      = dplyr::first(Effect_simple),
      label              = dplyr::first(label),
      .groups = "drop"
    )
}

if (!is.null(vars_genomic_plot) && !is.null(vars_exomic_plot)) {
  overlap_after <- sum(vars_genomic_plot$var_key %in% vars_exomic_plot$var_key)
}

##Case-Varianten Positionen (rote Linien im Density-Panel)
case_x <- c()
if (!is.null(vars_genomic_plot) && nrow(vars_genomic_plot) > 0) case_x <- c(case_x, vars_genomic_plot$scaled_x)
if (!is.null(vars_exomic_plot)  && nrow(vars_exomic_plot)  > 0) case_x <- c(case_x, vars_exomic_plot$scaled_x)
case_x <- case_x[is.finite(case_x)]


##Domänen mit Exons verknüpfen (für Plot)
##damit die funktionellen ELemente in den Exons angezeigt werden
if (!is.null(domain_df) && nrow(domain_df) > 0) {
  chr_name <- cds_df$chromosome_name[1]
  
  gr_exons <- GRanges(
    seqnames = chr_name,
    ranges   = IRanges(start = exon_df$exon_start,
                       end   = exon_df$exon_end),
    exon_id  = exon_df$exon_id
  )
  
  gr_domains <- GRanges(
    seqnames = chr_name,
    ranges   = IRanges(start = domain_df$genomic_start,
                       end   = domain_df$genomic_end),
    domain       = domain_df$domain,
    color        = domain_df$color,
    scaled_start = domain_df$scaled_start,
    scaled_end   = domain_df$scaled_end
  )
  
  hits <- findOverlaps(gr_domains, gr_exons)
  
  domain_exon_df <- data.frame(
    domain       = mcols(gr_domains)$domain[queryHits(hits)],
    color        = mcols(gr_domains)$color[queryHits(hits)],
    scaled_start = mcols(gr_domains)$scaled_start[queryHits(hits)],
    scaled_end   = mcols(gr_domains)$scaled_end[queryHits(hits)],
    exon_id      = mcols(gr_exons)$exon_id[subjectHits(hits)]
  )
  
  domain_df <- domain_exon_df
}

##################################################################################
###funktionelle Elemente

##funktionelle Elemente einlesen und integrieren
functional_elements_df <- NULL

if (!is.null(vars_genomic_noncoding) && nrow(vars_genomic_noncoding) > 0) {
  chr_name <- cds_df$chromosome_name[1]
  nc_gr <- GRanges(
    seqnames = chr_name,
    ranges = IRanges(
      start = vars_genomic_noncoding$pos,
      end   = vars_genomic_noncoding$pos
    )
  )
  
  functional_list <- list()
  
##ENCODE
if (file.exists(encode_file)) {
  encode_gr <- rtracklayer::import(encode_file)
  seqlevels(encode_gr) <- gsub("^chr", "", seqlevels(encode_gr))
  encode_gr <- encode_gr[seqnames(encode_gr) == chr_name]
  encode_hits <- subsetByOverlaps(encode_gr, nc_gr, ignore.strand = TRUE)
  
  if (length(encode_hits) > 0) {
    meta_enc <- as.data.frame(mcols(encode_hits))
    
    element_name <- if (ncol(meta_enc) >= 4) {
      as.character(meta_enc[[4]])
    } else if ("name" %in% names(meta_enc)) {
      as.character(meta_enc$name)
    } else {
      paste0("ENCODE_", seq_along(encode_hits))
    }
      
    encode_df <- data.frame(
      genomic_start = start(encode_hits),
      genomic_end   = end(encode_hits),
      element_name  = element_name,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::distinct(genomic_start, genomic_end, element_name, .keep_all = TRUE)
    
    encode_df$source <- "ENCODE"
    
    functional_list[["ENCODE"]] <- encode_df
  }
}
  
##RefSeq
if (file.exists(refseq_file)) {
  refseq_gr <- rtracklayer::import(refseq_file)
  seqlevels(refseq_gr) <- gsub("^chr", "", seqlevels(refseq_gr))
  refseq_gr <- refseq_gr[seqnames(refseq_gr) == chr_name]
  refseq_hits <- subsetByOverlaps(refseq_gr, nc_gr, ignore.strand = TRUE)
  
  if (length(refseq_hits) > 0) {
    meta_ref <- as.data.frame(mcols(refseq_hits))
    
    element_name <- if (ncol(meta_ref) >= 4) {
      as.character(meta_ref[[4]])
    } else if ("name" %in% names(meta_ref)) {
      as.character(meta_ref$name)
    } else {
      paste0("RefSeq_", seq_along(refseq_hits))
    }
      
    refseq_df <- data.frame(
      genomic_start = start(refseq_hits),
      genomic_end   = end(refseq_hits),
      element_name  = element_name,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::distinct(genomic_start, genomic_end, element_name, .keep_all = TRUE)
    
      refseq_df$source       <- "RefSeq"
    
    functional_list[["RefSeq"]] <- refseq_df
  }
}
  
  if (length(functional_list) > 0) {
    functional_elements_df <- dplyr::bind_rows(functional_list)
  }
}

if (!is.null(functional_elements_df) && nrow(functional_elements_df) > 0) {
  functional_elements_df <- functional_elements_df %>%
    mutate(
      scaled_start = map_genomic_to_scaled(genomic_start),
      scaled_end   = map_genomic_to_scaled(genomic_end),
      scaled_end   = ifelse(scaled_end == scaled_start,
                            scaled_end + 10,   
                            scaled_end)
    )
}

##################################################################################
###Erstellung des Plots

##Kombinierter Plot (Genom + Exom + Domänen)
shape_map <- c(
  Missense        = 16,  # gefüllter Kreis
  LOF             = 17,  # gefülltes Dreieck
  `Other/Unknown` = 15,  # gefülltes Quadrat
  `Non-coding`    = 18   # gefüllte Raute
  
)

shape_symbols <- c(
  LOF             = "▲",
  Missense        = "●",
  `Non-coding`    = "◆",
  `Other/Unknown` = "■"
)

shape_legend_text <- paste(
  paste(shape_symbols[names(shape_map)], names(shape_map)),
  collapse = ", "
)

##Gen-/Exom-Track: Exons, Domänen, coding Varianten
plot_combined_track <- function(vars_genomic_df,
                                vars_exomic_df,
                                vars_genomic_nc_plot, 
                                transcript_id,
                                strand,
                                domain_df = NULL,
                                phenotype_colors = NULL,
                                functional_elements_df = NULL) {
tx_min <- TX_MIN
tx_max <- TX_MAX

##dynamische Y-Achsenbeschriftung
has_nc <- !is.null(vars_genomic_nc_plot) && nrow(vars_genomic_nc_plot) > 0

if (has_nc) {
  y_breaks <- c(Y_GENOMIC_NC, Y_GENOMIC_LOF, Y_GENOMIC_MISS, Y_GENOMIC_OTHER,
                Y_EXOM_OTHER, Y_EXOM_MISS, Y_EXOM_LOF)
  y_labels <- c("non-coding (GS)", "LoF (GS)", "Missense (GS)", "Other (GS)",
                "Other (ES)", "Missense (ES)", "LoF (ES)")
} else {
  y_breaks <- c(Y_GENOMIC_LOF, Y_GENOMIC_MISS, Y_GENOMIC_OTHER,
                Y_EXOM_OTHER, Y_EXOM_MISS, Y_EXOM_LOF)
  y_labels <- c("LoF (GS)", "Missense (GS)", "Other (GS)",
                "Other (ES)", "Missense (ES)", "LoF (ES)")
}

##Genom- & Exom-Varianten aufbereiten
gen_plot <- NULL
snv_plot <- NULL
sv_plot  <- NULL

if (!is.null(vars_genomic_df) && nrow(vars_genomic_df) > 0) {
  gen_plot <- vars_genomic_df %>%
    dplyr::filter(!is.na(scaled_x)) %>%
    mutate(
      y_pos = dplyr::case_when(
        Effect_simple == "LOF"        ~ Y_GENOMIC_LOF,
        Effect_simple == "Missense"   ~ Y_GENOMIC_MISS,
        Effect_simple == "Non-coding" ~ Y_GENOMIC_NC,
        TRUE                          ~ Y_GENOMIC_OTHER
      ),
      x_cluster = floor(scaled_x / X_STACK_TOL)
    ) %>%
    group_by(x_cluster, y_pos) %>%
    arrange(scaled_x, .by_group = TRUE) %>%
    mutate(
      x_pos       = scaled_x,
      stack_index = dplyr::row_number(),
      y_pos       = y_pos + (stack_index - 1) * STACK_DY * 0.6
    ) %>%
    ungroup()
  
  sv_plot  <- gen_plot %>% dplyr::filter(is_sv)
  snv_plot <- gen_plot %>% dplyr::filter(!is_sv | is.na(is_sv))
}

ex_plot <- NULL
if (!is.null(vars_exomic_df) && nrow(vars_exomic_df) > 0) {
  ex_plot <- vars_exomic_df %>%
    dplyr::filter(!is.na(scaled_x)) %>%
    mutate(
      y_pos = dplyr::case_when(
        Effect_simple == "LOF"      ~ Y_EXOM_LOF,
        Effect_simple == "Missense" ~ Y_EXOM_MISS,
        TRUE                        ~ Y_EXOM_OTHER
      ),
      x_cluster = floor(scaled_x / X_STACK_TOL)
    ) %>%
    group_by(x_cluster, y_pos) %>%
    arrange(scaled_x, .by_group = TRUE) %>%
    mutate(
      x_pos       = scaled_x,
      stack_index = dplyr::row_number(),
      y_pos       = y_pos - (stack_index - 1) * STACK_DY * 0.6
    ) %>%
    ungroup()
}

##Plot aufbauen
p <- ggplot() 

separator_y <- c(
  if (has_nc) Y_GENOMIC_NC - 0.05,
  Y_GENOMIC_LOF   - 0.05,
  Y_GENOMIC_MISS  - 0.05,
  Y_GENOMIC_OTHER - 0.05,
  0               - 0.05,
  Y_EXOM_OTHER    - 0.05,
  Y_EXOM_MISS     - 0.05
)

p <- p +
  geom_segment(
    data = data.frame(y = separator_y),
    aes(
      x    = tx_min - INTRON_GAP,
      xend = tx_max + INTRON_GAP,
      y    = y,
      yend = y
    ),
    color        = "grey90",
    linewidth    = 0.4,
    inherit.aes  = FALSE
  ) +
  ## Gen-Basislinie bei y = 0
  geom_segment(
    aes(
      x    = tx_min - INTRON_GAP,
      xend = tx_max + INTRON_GAP,
      y    = 0,
      yend = 0
    ),
    linewidth   = GENE_LINE_WIDTH,
    inherit.aes = FALSE
  )

if (!is.null(snv_plot) && nrow(snv_plot) > 0) {
  p <- p +
    geom_segment(
      data = snv_plot,
      aes(x = x_pos, xend = x_pos,
          y = 0, yend = TICK_LEN),
      linewidth = 0.7,
      inherit.aes = FALSE
    ) 
}

if (!is.null(ex_plot) && nrow(ex_plot) > 0) {
  p <- p +
    geom_segment(
      data = ex_plot,
      aes(x = x_pos, xend = x_pos,
          y = 0, yend = -TICK_LEN),   
      linewidth = 0.7,
      inherit.aes = FALSE
    ) 
}

p <- p +
  geom_rect(
    data = exon_plot_df,
    aes(xmin = scaled_start, xmax = scaled_end,
        ymin = -EXON_HALF_HEIGHT, ymax = EXON_HALF_HEIGHT),
    fill = "grey97", 
    color = "black", 
    linewidth = EXON_BORDER_WIDTH,
    inherit.aes = FALSE) +
  
if (!is.null(vars_genomic_nc_plot) && nrow(vars_genomic_nc_plot) > 0) {
  
  func_ticks <- vars_genomic_nc_plot %>%
    mutate(x_pos = scaled_x)
  
  p <- p +
    geom_segment(
      data = func_ticks,
      aes(
        x    = x_pos,
        xend = x_pos,
        y    = Y_FUNC_RECT_TOP,
        yend = Y_FUNC_RECT_TOP + 0.05
      ),
      linewidth = 0.7,
      color     = "black",
      inherit.aes = FALSE
    )
}

if (!is.null(functional_elements_df) && nrow(functional_elements_df) > 0) {
  
  p <- p +
     geom_rect(
      data = functional_elements_df,
      aes(
        xmin = scaled_start,
        xmax = scaled_end,
        fill = source
      ),
      ymin        = Y_FUNC_RECT_BOTTOM,
      ymax        = Y_FUNC_RECT_TOP,
      alpha       = 0.6,
      color       = NA,
      inherit.aes = FALSE
    )
}

  if (!is.null(domain_df) && nrow(domain_df) > 0) {
    p <- p +
      geom_rect(
        data = domain_df,
        aes(xmin = scaled_start, xmax = scaled_end,
            ymin = -DOMAIN_HALF_HEIGHT, ymax = DOMAIN_HALF_HEIGHT,
            fill = domain),
        alpha = 0.6,
        color = NA,
        inherit.aes = FALSE
      )
  }  
    
 p <- p +
  geom_text(
    data = exon_plot_df,
    aes(x = (scaled_start + scaled_end)/2, y = 0, label = exon_id),
    size = EXON_LABEL_SIZE, 
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  )

##Gemeinsame Fill-Skala für Domänen + funktionelle Elemente
 
fill_values <- c()

#Domain-Farben aus domain_input übernehmen
if (!is.null(domain_input) && nrow(domain_input) > 0) {
 domain_colors <- setNames(domain_input$color, domain_input$domain)
 fill_values <- c(fill_values, domain_colors)
}

#Funktionelle Elemente automatisch erkennen
func_levels <- character(0)
func_labels <- character(0)

if (!is.null(functional_elements_df) && nrow(functional_elements_df) > 0) {
 func_levels <- unique(functional_elements_df$source)
 func_labels <- paste0(func_levels, "-Elements")
 
#Farben für functional elements
func_colors <- setNames(rep("#FF0000", length(func_levels)), func_levels)
fill_values <- c(fill_values, func_colors)
}
 
#Reihenfolge/Labels automatisch bauen
domain_levels <- if (!is.null(domain_input) && nrow(domain_input) > 0) domain_input$domain else character(0)
domain_labels <- gsub("_", " ", domain_levels)

all_levels <- c(domain_levels, func_levels)
all_labels <- c(domain_labels, func_labels)
 
#Scale nur setzen, wenn etwas da ist
present_levels <- intersect(all_levels, names(fill_values))

if (length(present_levels) > 0) {
 p <- p + scale_fill_manual(
   values = fill_values[present_levels],
   name   = "Additional elements",
   breaks = present_levels,
   labels = all_labels[match(present_levels, all_levels)]
 )
}
 
##GENOM-Varianten (oben)
##SNV-Plot
  if (!is.null(snv_plot) && nrow(snv_plot) > 0) {
    p <- p +
      geom_point(
        data = snv_plot,
        aes(x = x_pos, y = y_pos,
            color = Phenotype_complete,
            shape = Effect_simple,
            size  = in_both)
      )
    
    snv_label <- snv_plot %>%
      dplyr::filter(!is.na(label), !is.na(x_pos), !is.na(y_pos))
  
  if (nrow(snv_label) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = snv_label,
          aes(x = x_pos, y = y_pos, label = label),
          angle         = LABEL_ANGLE,
          nudge_x       = LABEL_NUDGE_X,
          vjust         = 0,
          hjust         = 0.5,
          segment.color = NA,
          max.overlaps  = 50,
          size          = 5
        )
    }
  }
  
##SV-Plot 
  if (nrow(sv_plot) > 0) {
    p <- p +
      geom_rect(
        data = sv_plot,
        aes(
          xmin  = scaled_start,
          xmax  = scaled_end,
          color = Phenotype_complete,
         ),
        fill  = NA,
        ymin  = Y_GENOMIC_OTHER - SV_HALF_HEIGHT,
        ymax  = Y_GENOMIC_OTHER + SV_HALF_HEIGHT,
        linewidth = 1.0,
        show.legend = FALSE
      )
      
    sv_label <- sv_plot %>%
      mutate(
        label_x = (scaled_start + scaled_end) / 2,
        label_y = Y_GENOMIC_OTHER + 0.04
      ) %>%
      dplyr::filter(!is.na(label))
      
    if (nrow(sv_label) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = sv_label,
          aes(x = label_x, y = label_y, label = label),
          angle         = LABEL_ANGLE,
          nudge_x       = LABEL_NUDGE_X,
          vjust         = 0,
          hjust         = 0.5,
          segment.color = NA,
          max.overlaps  = 50,
          size          = 5
        )
    }
  }
  
##EXOM-Varianten (unten)
##Exom-Plot
if (!is.null(ex_plot) && nrow(ex_plot) > 0) {
  p <- p +
    geom_point(
      data = ex_plot,
      aes(x = x_pos, y = y_pos,
          color = Phenotype_complete,
          shape = Effect_simple),
      size = 3
    )
  
  ex_label <- ex_plot %>%
    dplyr::filter(!is.na(label), !is.na(x_pos), !is.na(y_pos))

    if (nrow(ex_label) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = ex_label,
          aes(x = x_pos, y = y_pos, label = label),
          angle         = LABEL_ANGLE,
          nudge_x       = LABEL_NUDGE_X, 
          vjust         = 0,
          hjust         = 0.5,
          segment.color = NA,
          max.overlaps  = 50,
          size          = 5
        )
    }
  }
  
##Skalen, Legenden & Achsen
p <- p + scale_shape_manual(values = shape_map, guide = "none")

p <- p + scale_size_manual(
  values = c(`FALSE` = 3, `TRUE` = 4),  # 3 = nur Genom, 5 = Genom+Exom
  guide  = "none" 
)
  
if (!is.null(phenotype_colors)) {
  p <- p + scale_color_manual(values = phenotype_colors,
                              na.value = "black",
                              name = "")
} else {
  p <- p + scale_color_discrete(name = "")
}

p <- p + ggtitle(plot_title)       

p <- p +
  theme_cowplot() +
  theme(
    legend.position      = "bottom",
    legend.justification = "center",
    legend.box           = "vertical",
    
    legend.title = element_text(size = 16),   # ← HIER
    legend.text  = element_text(size = 16),   # hast du schon
    
    axis.line.x  = element_line(color = "black", linewidth = 0.6),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y  = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y  = element_text(size = 16, hjust = 1),
    axis.title.y = element_blank(),
    plot.title   = element_text(size = 20, face = "bold")
  )+

  xlab(paste0("Transcript: ", transcript_id,
              "    Shapes: ", shape_legend_text)) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    limits = c(-0.50, 0.65)
  ) +
  guides(
    color = guide_legend(override.aes = list(
      shape = 16,
      size  = 2,
      fill  = NA
    ))
  )

  if (strand == -1) {
    p <- p + scale_x_reverse(
      limits = c(tx_max + INTRON_GAP, tx_min - INTRON_GAP)
    )
  } else {
    p <- p + scale_x_continuous(
      limits = c(tx_min - INTRON_GAP, tx_max + INTRON_GAP)
    )
  }

p
}

##Density-Plot UKB
plot_density_panel <- function(x_positions, tx_min, tx_max,
                               strand = 1, vlines = NULL,
                               n = 1024, adjust = 1) {
  
  x_positions <- as.numeric(unlist(x_positions, use.names = FALSE))
  x_positions <- x_positions[is.finite(x_positions)]
  
  if (length(x_positions) < 2) {
    return(ggplot() + theme_void())
  }
  
  dens <- density(
    x_positions,
    from = tx_min - INTRON_GAP,
    to   = tx_max + INTRON_GAP,
    n = n,
    adjust = adjust
  )
  
  density_df <- data.frame(x = dens$x, y = dens$y)
  
  p <- ggplot(density_df, aes(x = x, y = y)) +
    geom_line(linewidth = 0.7) +
    labs(y = "UK_Biobank\nrelative density\n(AF<0.1%)", x = NULL) +
    theme_cowplot() +
    theme(
      plot.title   = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  if (!is.null(vlines) && length(vlines) > 0) {
    vdf <- data.frame(x = vlines[is.finite(vlines)])
    p <- p + geom_vline(
      data = vdf, 
      aes(xintercept = x), 
      linewidth = 0.6, 
      color = "red",
      alpha = 0.8)
  }
  
  if (strand == -1) {
    p <- p + scale_x_reverse(limits = c(tx_max + INTRON_GAP, tx_min - INTRON_GAP))
  } else {
    p <- p + scale_x_continuous(limits = c(tx_min - INTRON_GAP, tx_max + INTRON_GAP))
  }
  
  p
}

##Density-Plot gnomAD
 plot_density_panel2 <- function(x_positions2, tx_min, tx_max,
                               strand = 1, vlines = NULL,
                               n = 1024, adjust = 1) {
  
  x_positions2 <- as.numeric(unlist(x_positions2, use.names = FALSE))
  x_positions2 <- x_positions2[is.finite(x_positions2)]
  
  if (length(x_positions2) < 2) {
    return(ggplot() + theme_void())
  }
  
  dens2 <- density(
    x_positions2,
    from = tx_min - INTRON_GAP,
    to   = tx_max + INTRON_GAP,
    n = n,
    adjust = adjust
  )
  
  density_df2 <- data.frame(x = dens2$x, y = dens2$y)
  
  p <- ggplot(density_df2, aes(x = x, y = y)) +
    geom_line(linewidth = 0.7) +
    labs(y = "gnomAD 4.0\nrelative density\n(AF<0.1%)", x = NULL) +
    theme_cowplot() +
    theme(
      plot.title   = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  if (!is.null(vlines) && length(vlines) > 0) {
    vdf2 <- data.frame(x = vlines[is.finite(vlines)])
    p <- p + geom_vline(
     data = vdf2, 
      aes(xintercept = x), 
      linewidth = 0.6, 
      color = "red",
      alpha = 0.8)
     }
  
  if (strand == -1) {
   p <- p + scale_x_reverse(limits = c(tx_max + INTRON_GAP, tx_min - INTRON_GAP))
 } else {
   p <- p + scale_x_continuous(limits = c(tx_min - INTRON_GAP, tx_max + INTRON_GAP))
 }
  
 p
}


##Fehlende Phänotyp-Farben automatisch ergänzen
if (!is.null(phenotype_colors)) {
  all_pheno <- c()
  if (!is.null(vars_genomic_plot))    all_pheno <- c(all_pheno, vars_genomic_plot$Phenotype_complete)
  if (!is.null(vars_exomic_plot))     all_pheno <- c(all_pheno, vars_exomic_plot$Phenotype_complete)
  if (!is.null(vars_genomic_nc_plot)) all_pheno <- c(all_pheno, vars_genomic_nc_plot$Phenotype_complete)
  
  all_pheno <- unique(all_pheno[!is.na(all_pheno)])
  
  missing_pheno <- setdiff(all_pheno, names(phenotype_colors))
  
  if (length(missing_pheno) > 0) {
    extra_cols <- grDevices::rainbow(length(missing_pheno))
    names(extra_cols) <- missing_pheno
    phenotype_colors <- c(phenotype_colors, extra_cols)
  }
}

##################################################################################
###Plot ausführen

##Plot zeichnen
if ((is.null(vars_genomic_plot)    || nrow(vars_genomic_plot)    == 0) &&
    (is.null(vars_exomic_plot)     || nrow(vars_exomic_plot)     == 0) &&
    (is.null(vars_genomic_nc_plot) || nrow(vars_genomic_nc_plot) == 0)) {
  stop("Es wurden keine Varianten gefunden – nichts zu plotten.")
}

##Gen-Plot (nur coding-Varianten)
gene_plot <- plot_combined_track(
  vars_genomic_df        = vars_genomic_plot,
  vars_exomic_df         = vars_exomic_plot,
  vars_genomic_nc_plot   = vars_genomic_nc_plot,
  transcript_id          = transcript_id,
  strand                 = strand,
  domain_df              = domain_df,
  functional_elements_df = functional_elements_df,
  phenotype_colors       = phenotype_colors
)

## Gemeinsames y-Maximum für beide Density-Panels bestimmen
dens_ukb  <- if (length(UKB_scaled_x) >= 2) {
  density(UKB_scaled_x,
          from = TX_MIN - INTRON_GAP, to = TX_MAX + INTRON_GAP,
          n = DENSITY_N, adjust = DENSITY_ADJUST)
} else NULL

dens_gnom <- if (length(gnomAD_scaled_x) >= 2) {
  density(gnomAD_scaled_x,
          from = TX_MIN - INTRON_GAP, to = TX_MAX + INTRON_GAP,
          n = DENSITY_N, adjust = DENSITY_ADJUST)
} else NULL

YMAX_DENS <- max(
  if (!is.null(dens_ukb))  dens_ukb$y  else 0,
  if (!is.null(dens_gnom)) dens_gnom$y else 0,
  na.rm = TRUE
)

density_plot <- plot_density_panel(
  x_positions = UKB_scaled_x,     
  tx_min = TX_MIN,
  tx_max = TX_MAX,
  strand = strand,
  vlines = case_x,                
  n = DENSITY_N,
  adjust = DENSITY_ADJUST
) + coord_cartesian(ylim = c(0, YMAX_DENS))

density_plot2 <- plot_density_panel2(
  x_positions2 = gnomAD_scaled_x,     
  tx_min = TX_MIN,
  tx_max = TX_MAX,
  strand = strand,
  vlines = case_x,                
  n = DENSITY_N,
  adjust = DENSITY_ADJUST
) + coord_cartesian(ylim = c(0, YMAX_DENS))

final_plot <- cowplot::plot_grid(
  gene_plot,
  ggplot() + theme_void(),
  density_plot,
  density_plot2,
  ncol = 1,
  align = "v",
  rel_heights = c(3.6, 0.5, 0.6, 0.6)
)


print(final_plot)

##################################################################################
###Plot speichern
ggsave(
  filename = paste0(gene_name, "_genomic_plot.png"), 
  plot     = final_plot,
  width    = 35,
  height   = 15,
  dpi      = 300,
  device   = ragg::agg_png,
  bg       = "white"
)

