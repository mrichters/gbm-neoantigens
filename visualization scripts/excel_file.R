# 11.5.19

# import libraries
library(openxlsx)
library(readxl)

# set wd
input_dir <- "~/Desktop/GBM/immune_infiltration/danaher_analysis"
setwd(input_dir)

# danaher genes file
danaher_genes <- paste(input_dir, "danaher_gene_set.xlsx", sep='/')

# list all expression files in directory
expression_files <- list.files(paste(input_dir, "danaher_files", sep='/'))
# make the files into a list
#expression_data <- lapply(expression_files, read.table, header = TRUE, sep = '\t')
# initialize final df
final_df <- data.frame(expression_data[[1]][, c(1)])
# initialize vector for column names of final_df
#column_names <- c(names(final_df))
column_names <- c()
# for loop over expression_files
for ( i in expression_files ) {
  # read tsv file
  data <- read.table(i, header = TRUE, sep = '\t')
  # get abundance info
  tpm_column <- data[, 3]
  # add abundance info to final_df
  final_df <- cbind(final_df, tpm_column)
  # create column name
  header <- unlist(strsplit(i, ".gene"))[1]
  # add header to column names vector
  column_names <- c(column_names, header)
} 
# add column names to final_df
names(final_df) <- c("Cell Marker", "GeneSymbol", column_names)
# save final_df to excel file
#WriteXLS(final_df, "transcript_expression.xlsx", SheetNames = "all_transcripts")
# cibersort input file
write.table(final_df, file = "GBM_estimate_mixture_file.txt", quote = FALSE, sep = '\t', row.names = FALSE)
