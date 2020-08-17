# 09.26.19

library(readxl)
library(dplyr)
library(rsed)

df <- read.csv("~/Desktop/GBM/junction_pvalues_default.tsv", sep = '\t') 
df.subset <- df %>% filter(anchor != "N" & anchor != "DA" & pvalue <= 0.05)
df.subset.sort <- df.subset[order(df.subset$pvalue), ]

samples_vector <- c()
for ( sample_list in df.subset.sort$samples ) {
  if ( grepl("c", sample_list) ) {
    new_name <- substr(sample_list, 3, as.numeric(nchar(sample_list))-1)
    samples_vector <- c(samples_vector, new_name)
  } else {
    samples_vector <- c(samples_vector, sample_list)
  }
}
df.subset.sort$samples <- samples_vector

write.table(df.subset.sort, "~/Desktop/GBM/junction_pvalues_default_priority.tsv", sep = "\t", quote=F, row.names = F)
