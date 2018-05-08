#!/usr/bin/env nextflow

cancers = ["OV", "BRCA"]

/*
    Download TCGA data, perform regression and save results as csv
*/

process fit_cnv_expression {
    cache true

    publishDir "$baseDir/data/gene_tables", mode: 'copy',
        saveAs: {filename -> "tcga_gene_table_${cancer}.csv"}

    input:
    each cancer from cancers

    output:
    set val(cancer), file("tcga_gene_table.csv") into tcga_gene_tables

    """
    Rscript $baseDir/scripts/tcga_gene_selection.R \
        --cancer $cancer
    """
}

process graph_cnv_expression {
    cache false
    
    publishDir "$baseDir/figs/${cancer}", mode: 'copy'

    input:
    set val(cancer), file("tcga_gene_table.csv") from tcga_gene_tables

    output:
    file("*.png")

    """
    Rscript $baseDir/scripts/graph_tcga_gene_selection.R \
    --cancer $cancer \
    --input_csv tcga_gene_table.csv
    """
}