library(biomaRt)
library(tidyverse)
library(chunkR)
library(readxl)
library(xlsx)

setwd("~/Desktop/GBM/cnvkit/")

# ensembl_99 <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
# file_chunks <- chunker("unique_total_predicted_genes.tsv", chunksize = 500, has_rownames = F)
# 
# write_lines("ensembl_gene_id\texternal_gene_name", "converted_gene_ids.tsv", append = FALSE)
# 
# repeat {
#     if ( next_chunk(file_chunks) == F ) {
#         break
#     }
#     df <- get_table(file_chunks)
#     gene_df <- getBM(c("ensembl_gene_id", "external_gene_name"), filters = 'ensembl_gene_id', values = df$ensembl_gene_ids, mart = ensembl_99)
#     for (row in 1:nrow(gene_df)) {
#         file_line <- paste(c(gene_df$ensembl_gene_id[row], gene_df$external_gene_name[row]), collapse='\t')
#         write_lines(file_line, "converted_gene_ids.tsv", append = TRUE)
#     }
#     print(get_completed(file_chunks))
# }

#### append gene names to gene ids' excel sheet ####

conversion_df <- read.table("converted_gene_ids.tsv", header=T)
conversion_df <- mutate_all(conversion_df, as.character)

sample_names <- excel_sheets("trusted_genes_0.4.xlsx")

wb = createWorkbook()

for ( sample in sample_names ) {
    excel_sheet <- as.data.frame(read_excel("trusted_genes_0.4.xlsx", sheet = sample))
    write_df <- c()
    for ( row in 1:nrow(excel_sheet) ) {
        gene_id <- excel_sheet[row, sample]
        if ( !dim(filter(conversion_df, ensembl_gene_id == gene_id))[1] == 0 ) {
            gene_name <- filter(conversion_df, ensembl_gene_id == gene_id)$external_gene_name
            write_line <- data.frame("Ensembl_Gene_ID" = gene_id, "Gene_Name" = gene_name)
            write_df <- rbind(write_df, write_line)
        }
    }
    sheet = createSheet(wb, sample)
    if (sample == sample_names[1] ) {
        addDataFrame(write_df, sheet=sheet, startColumn = 1, row.names = F, col.names = T)
    } else {
        addDataFrame(write_df, sheet=sheet, startColumn = 1, row.names = F, col.names = T, append = T)
    }
}

saveWorkbook(wb, "cnvkit_significant_genes_0.4.xlsx")

