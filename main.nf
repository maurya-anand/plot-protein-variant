#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    def gene_list_ch = channel.fromPath(params.gene_list)
        .splitText()
        .map { gene -> gene.trim() }

    def info_ch = FETCH_GENE_INFO(gene_list_ch)

    def genomic_ch = channel.fromPath(params.genomic)

    def exonic_ch = channel.fromPath(params.exonic)

    def gnomad_ch = channel.fromPath(params.gnomAD)

    def uk_biobank_ch = channel.fromPath(params.uk_Biobank)

    def colour_phenotypes_ch = channel.fromPath(params.colour_phenotypes)

    def genomic_noncoding_ch = params.genomic_noncoding
        ? channel.fromPath(params.genomic_noncoding)
        : channel.of(file('NO_FILE_GENOMIC_NONCODING'))

    def sv_ch = params.sv
        ? channel.fromPath(params.sv)
        : channel.of(file('NO_FILE_SV'))

    def encode_file_ch = params.encode_file
        ? channel.fromPath(params.encode_file)
        : channel.of(file('NO_FILE_ENCODE'))

    def refseq_file_ch = params.refseq_file
        ? channel.fromPath(params.refseq_file)
        : channel.of(file('NO_FILE_REFSEQ'))

    def plot_input_ch = info_ch.info
        .combine(genomic_ch)
        .combine(genomic_noncoding_ch)
        .combine(exonic_ch)
        .combine(sv_ch)
        .combine(encode_file_ch)
        .combine(refseq_file_ch)
        .combine(gnomad_ch)
        .combine(uk_biobank_ch)
        .combine(colour_phenotypes_ch)

    PLOT_VARIANTS(plot_input_ch)
}

process FETCH_GENE_INFO {
    tag "${id}"
    publishDir ("${params.outdir}/${id}"), mode: 'copy'
    publishDir ("${params.outdir}/${id}/FETCH_GENE_INFO_log"), mode: 'copy', pattern: ".command.log"

    input:
    val id

    output:
    tuple val(id), path("${id}_domain.tsv"), path("${id}_transcript.tsv"), emit: info
    path ".command.log", optional: true

    script:
    """
    gene-to-protein-domains.py --gene ${id}
    """
}

process PLOT_VARIANTS {
    tag "${id}"
    publishDir ("${params.outdir}/${id}"), mode: 'copy'
    publishDir ("${params.outdir}/${id}/PLOT_VARIANTS_log"), mode: 'copy', pattern: ".command.log"

    input:
    tuple val(id), path("${id}_domain.tsv"), path("${id}_transcript.tsv"), path(genomic_file), path(genomic_noncoding_file), path(exonic_file), path(sv_file), path(encode_file), path(refseq_file), path(gnomAD_file), path(UK_Biobank_file), path(colour_phenotypes_file)

    output:
    path ("${id}_genomic_plot.png"), emit: plot
    path ".command.log", optional: true

    script:
    """
    export HOME="\$PWD"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME/R/biomaRt"    
    rm NO_FILE*
    ls -ltrh
    transcript=\$(grep '^canonical_transcript' ${id}_transcript.tsv | cut -f2)
    Plots_Genomic_Exonic_Protein_Effect_sorted.R \
        --gene_name ${id} \
        --gene_domain ${id}_domain.tsv \
        --transcript_id \$transcript \
        --genomic ${genomic_file} \
        --genomic_noncoding ${genomic_noncoding_file} \
        --exonic ${exonic_file} \
        --sv ${sv_file} \
        --encode_file ${encode_file} \
        --refseq_file ${refseq_file} \
        --gnomAD ${gnomAD_file} \
        --UK_Biobank ${UK_Biobank_file} \
        --colour_phenotypes ${colour_phenotypes_file} \
    """
}
