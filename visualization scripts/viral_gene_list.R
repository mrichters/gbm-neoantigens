# 11.15.19
# identify HPV viral genes from file from gdc.cancer.gov

# import packages
library(tidyverse)
library(readxl)
library(openxlsx)

# set wd
setwd("~/Desktop/GBM/targeted_resequencing")

# read in file
df <- read_tsv("GRCh83.d1.vd1_virus_decoy.txt")
df$all_NA <- NA
for ( row in 1:nrow(df) ) {
  df[row, "all_NA"] <- all(is.na(df[row, ]))
}
df.filter <- filter(df, all_NA == FALSE)

total_gene_list <- as.data.frame(df.filter$Abbreviation)
colnames(total_gene_list) <- "Viral Genes List"
write_tsv(total_gene_list, "viral_genes_list.tsv")
