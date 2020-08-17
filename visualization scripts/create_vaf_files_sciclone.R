library(hash)
library(tidyverse)

setwd("~/Desktop/GBM/GBM065/")

tsv_files <- list.files("vtt_final/")
for ( file in tsv_files ) {
    # read in file
    print(file)
    readcounts <- read.table(paste("vtt_final/", file, sep = ""), header = TRUE, sep = '\t')
    # dictionaries with sample name replacements
    gsub_replacements <- hash(keys=c("TUMOR.variant", "TUMOR.variant2", "TUMOR.variant3", "TUMOR.variant4", "TUMOR.variant5", "TUMOR.variant6"), values=c("GBM065-primary", "GBM065-re1", "GBM065-re2", "GBM065-re3", "GBM065-re4", "GBM065-re5"))
    if ( grepl("primary", file) ) {
        name <- "GBM065-primary"
        readcounts$Sample <- name
        colnames(readcounts) <- gsub("TUMOR.variant", name, colnames(readcounts))
        colnames(readcounts) <- gsub(paste(name, "\\.", sep=""), "", colnames(readcounts))
    } else {
        hash_key <- unlist(str_split(file, "_"))[2]
        readcounts$Sample <- gsub_replacements[[ hash_key ]]
        colnames(readcounts) <- gsub(paste(hash_key, "\\.", sep=""), "", colnames(readcounts))
    }
    readcounts$AD <- as.character(readcounts$AD)
    # add alt reads and VAF for final file
    readcounts$`NAD` <- sapply(readcounts$AD, function(x) return(unlist(strsplit(x, ','))[1]))
    readcounts$`TAD` <- sapply(readcounts$AD, function(x) return(unlist(strsplit(x, ','))[2]))
    readcounts$`NAD` <- as.numeric(readcounts$`NAD`)
    readcounts$`TAD` <- as.numeric(readcounts$`TAD`)
    readcounts <- transform(readcounts, VAF = (TAD / DP)*100)
    readcounts$VAF[is.nan(readcounts$VAF)] <- 0
    # now ready for final files: chr	pos	ref_reads	var_reads	vaf
    readcounts.final <- readcounts[c("CHROM", "POS", "NAD", "TAD", "VAF")]
    colnames(readcounts.final) <- c("chr", "pos", "ref_reads", "var_reads", "vaf")
    # write to file
    new_name <- paste(unique(readcounts$Sample), ".tsv", sep="")
    print(new_name)
    write_tsv(readcounts.final, new_name)
}
# these files go in sciclone_vaf_files/
