# import packages
library(tidyverse)

# set wd
setwd("~/Desktop/GBM")

unique_ensg <- read_tsv("unique_total_predicted_genes.tsv")

gene_files <- list.files("trusted_genes_files/")
#df_list <- lapply(gene_files, function(a) sub(".trusted-genes.txt", "", a))


#mylist <- vector("list", 59)
#i <- 1
for ( file in gene_files ) {
  name <- sub(".trusted-genes.txt", "", file)
  print(name)
  single_file <- as.data.frame(read_tsv(paste("trusted_genes_files/", file, sep=""), col_names = F))
  colnames(single_file) <- "ENSEMBL Gene IDs"
  ensg_row <- data.frame(name = single_file$`ENSEMBL Gene IDs`)
  #mylist[[i]] <- ensg_row
  #i <- i + 1
}

