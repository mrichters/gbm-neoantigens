setwd("~/Desktop/GBM/htseq-count/")
tsvs <- list.files("~/Desktop/GBM/htseq-count/")

df <- c()
for (f in tsvs) {
    sample <- gsub(".gene_counts.tsv", "", f)
    counts <- as.data.frame(read_tsv(f, col_names = c("Gene", sample)))
    if ( is.null(dim(df)) ) {
        df <- rbind(df, counts)
    } else {
        merged <- merge(df, counts, by="Gene")
        df <- merged
    }
}

library(WriteXLS)
WriteXLS(df, "htseq-counts_final.xlsx")
