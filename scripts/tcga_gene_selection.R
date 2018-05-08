
library(aargh)
library(TCGA2STAT)
library(tidyverse)
library(broom)


tcga_gene_selection <- function(cancer = "OV") {
  
  rnaseq.ov <- getTCGA(disease = cancer, data.type = "RNASeq", type = "RPKM", clinical = TRUE)
  cnvsnp <- getTCGA(disease = cancer, data.type = "CNV_SNP")
  
  merged_data <- OMICSBind(rnaseq.ov$dat, cnvsnp$dat)
  
  rna_seq <- merged_data$X
  cnv <- merged_data$Y
  
  rownames(rna_seq) <- sub("d1.", "", rownames(rna_seq), fixed = TRUE)
  rownames(cnv) <- sub("d2.", "", rownames(cnv), fixed = TRUE)
  
  common_genes <- intersect(rownames(rna_seq), rownames(cnv))
  
  rna_seq <- rna_seq[common_genes,]
  cnv <- cnv[common_genes,]
  cnv <- as.matrix(cnv)
  
  genes <- rownames(rna_seq)
  
  fit_tidy_model <- function(gene) {
    y <- log(rna_seq[gene,] + 1)
    x <- cnv[gene,]
    
    fit <- lm(y ~ x)
    glance(fit) %>% mutate(effect_size = coef(fit)[2])
  }
  
  df <- map_df(genes, fit_tidy_model)
  df <- as_tibble(df)
  df <- mutate(df, gene = genes, q.value = p.adjust(p.value))
  
  df <- select(df, gene, everything()) %>% 
    arrange(q.value)

  write_csv(df, "tcga_gene_table.csv")
  
}

aargh(tcga_gene_selection)