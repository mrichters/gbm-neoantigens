# 09.09.19
## CREATE IMMUNE HEATMAP ##

# import packages
library(readxl)
library(tidyverse)
library(ggsci)
library(gridExtra)
library(gtable)
library(scales)
library(ComplexHeatmap)
library(textshape)
library(circlize)

# set wd
setwd("~/Desktop/GBM/immune_infiltration/")

# identify genes
#immune_df <- as.data.frame(read_excel("~/Desktop/GBM/ncomms3612-s2.xlsx", sheet = "gene"))
#gene_list <- filter(immune_df, Set == "Immune141_UP")$Gene

# determine RNA expression per sample per gene
# file I need: kallisto gene output
# create samples and tumors list
samples <- excel_sheets(path="immune_genes_expression.xlsx")
tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[1])))

total_df <- data.frame()
for ( tumor in tumors ) {
  tumor_df <- c()
  for ( name in samples[grep(tumor, samples)] ) {
    expression <- as.data.frame(read_excel("immune_genes_expression.xlsx", sheet = name))
    categories <- unique(expression$category)
    for ( item in categories ) {
      lines <- filter(expression, category == item)
      genes_list <- unique(lines$gene_name)
      for ( single_gene in genes_list ) {
        gene_lines <- filter(lines, gene_name == single_gene)
        if ( nrow(gene_lines) > 1 ) {
          abundance <- mean(gene_lines$abundance)
          average_df_line <- cbind(gene_lines[1, c(1:3)], abundance, gene_lines[1, c(6,7)])
          tumor_df <- rbind(tumor_df, average_df_line)
        } else {
          tumor_df <- rbind(tumor_df, gene_lines[, c(1,2,3,5,6,7)])
        }
      }
    }
    #expression_df <- filter(expression[,c(1,3)], gene_name %in% gene_list)
    #expression_df.uniq <- expression_df[!duplicated(expression_df[,"gene_name"]),]
    #expression_df.order <- expression_df.uniq[order(expression_df.uniq$gene_name), ]
    sample_df <- data.frame(tumor_df$sample, tumor_df$gene_name, log(tumor_df$abundance))#as.numeric(log2(tumor_df$abundance)))
    sample_df[sample_df == "-Inf"] <- 0
    colnames(sample_df) <- c("Sample", "Gene Name", "Abundance")
  }
  total_df <- rbind(total_df, sample_df)
}

total_df.format <- spread(total_df, key = Sample, value = Abundance)#key = sample, value = abundance))
rownames(total_df.format) <- total_df.format[,1]
total_df.format <- total_df.format[, -c(1)]
total_df.format <- total_df.format[c("CD8A", "CD8B", "CD4", "PDCD1", "CD274", "CD3D", "CD3E", "CD3G"), ]
total_matrix_brmet <- total_df.format[, c(1:30)]
total_matrix_gbm

total_matrix <- data.matrix(total_df.format)
col_fun = colorRamp2(c(-2, 0, 2), c("red", "white", "blue"))
col_fun(seq(-3, 3))
#Heatmap(total_matrix, row_order = c("CD4","CD8A", "CD8B","CD3D", "CD3E", "CD3G", "PDCD1", "CD274"), name = "log(TPM)", width = unit(31, "cm"), height = unit(12, "cm"))

#Heatmap(total_matrix, name = "log(TPM)", width = unit(31, "cm"), height = unit(25, "cm"), column_order = samples)
Heatmap(total_matrix, name = "log(TPM)", width = unit(31, "cm"), height = unit(15, "cm"), column_order = samples)
