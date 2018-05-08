library(tidyverse)
library(glue)
library(aargh)

theme_set(theme_bw(base_size = 11))

graph_tcga_gene_selection <- function(cancer = "OV", input_csv = "input.csv") {
  
  df <- read_csv(input_csv)
  
  volcano_plot <- ggplot(df, aes(x = effect_size, y = r.squared)) +
    geom_point(alpha = 0.1) +
    labs(x = "Effect size, lm(expression ~ cnv)",
         y = expression(R^2))
  
  qval_histo <- qplot(df$q.value, bins = 100) +
    geom_vline(xintercept = 0.05, color = 'red', linetype = 2) +
    labs(x = "FDR adjusted p value")
  
  ggsave(glue("R2_effect_{cancer}.png"), volcano_plot, width = 5, height = 4)
  
  ggsave(glue("qval_histogram_{cancer}.png"), qval_histo, width = 6, height = 3)
  
}

aargh(graph_tcga_gene_selection)