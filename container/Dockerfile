FROM rocker/r-ver:4

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    ca-certificates \
    fonts-dejavu \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled --ncpus -1 \
    dplyr readxl stringr ggplot2 cowplot ggrepel ragg readr BiocManager \
    optparse colorspace \
    && rm -rf /tmp/downloaded_packages

RUN R -q -e 'BiocManager::install(c("rtracklayer", "GenomicRanges", "biomaRt"))'