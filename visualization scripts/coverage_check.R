# 08.24.20
# merge exome and targeted resequencing variant readcounts into one file
# Rscript coverage_format_files.R <sample_name> <exome_file> <resequencing_file> <output_dir>

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

format_files <- function(sample, path, seq_type) {
    data <- read.table(path, header=T)
    data <- data %>% select(-contains("GT"))
    AD_col <- names(data)[grep("AD", names(data))]
    s <- separate(data, AD_col, c(paste0(AD_col, ".REF"), paste0(AD_col, ".ALT")), sep=',', convert=T)
    VAF_col <- gsub("AD", "VAF", AD_col); DP_col <- gsub("AD", "DP", AD_col)
    s[VAF_col] <- s[paste0(AD_col, ".ALT")] / s[DP_col]
    sample_names <- sapply(names(s)[5:length(names(s))], function(x) gsub(paste0(gsub("-", ".", sample), '.'), "", x))
    new_names <- sapply(sample_names, function(x) return(paste(x, seq_type, sep='.')))
    names(s) <- c(names(s)[1:4], new_names)
    s[is.na(s)] <- 0
    return(s)
}

# format df from file
exome_df <- format_files(args[1], args[2], "exome")
ts_df <- format_files(args[1], args[3], "ts")

# merge df
c <- merge(exome_df, ts_df, by=names(exome_df)[1:4])
c <- c %>% arrange(CHROM, POS)

# filter variant calls
c_fil <- c %>% filter((DP.exome >= 30 & DP.ts >= 50 & VAF.exome >= 0.025 & VAF.ts >= 0.025) | (DP.exome >= 30 & DP.ts < 50 & VAF.exome >= 0.025))

# combine exome and ts readcounts for all passing variants
c_fil$DP <- c_fil$DP.exome + c_fil$DP.ts
c_fil$AD.REF <- c_fil$AD.REF.exome + c_fil$AD.REF.ts
c_fil$AD.ALT <- c_fil$AD.ALT.exome + c_fil$AD.ALT.ts
c_fil$VAF <- c_fil$AD.ALT / c_fil$DP
# subset df
c_fil_final <- c_fil %>% select(-contains(c("exome", "ts")))

# write merged df to file
write.table(c_fil_final, paste0(args[4], args[1], "_filtered.annotated.tsv"), sep='\t', quote=F, row.names=F)