library(WriteXLS)

setwd("~/Desktop/GBM/tsv/")

# for multiple folders of results tsvs
for ( folder in list.files() ) {
    df <- data.frame()
    for ( tsv in list.files(folder) ) {
        tsv_df <- read.table(paste(folder, "/", tsv, sep=""), header=T, sep='\t')
        sample_name <- gsub(".RNA.*\\.[a-z].*\\.tsv", "", tsv)
        print(sample_name)
        df_add <- data.frame("samples" = sample_name, tsv_df)
        #head(df_add)
        df <- rbind(df, df_add)
    }
        df.spread <- spread(df, samples, tpm)
        WriteXLS(df.spread, paste("../", folder, ".xlsx", sep=""))
}

# for a single output folder
df <- data.frame()
for ( tsv in list.files() ) {
    tsv_df <- read.table(paste(tsv, sep=""), header=T, sep='\t')
    sample_name <- gsub(".RNA.*\\.[a-z].*\\.tsv", "", tsv)
    print(sample_name)
    df_avg <- aggregate(tsv_df$abundance, list(tsv_df$gene_name), mean)
    names(df_avg) <- c("gene_name", "abundance")
    df_avg$sample <- sample_name
    df <- rbind(df, df_avg)
}
df.spread <- spread(df, gene_name, abundance)
WriteXLS(df.spread, paste("../", folder, ".xlsx", sep=""))
