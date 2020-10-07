# 08.31.20
# creating a df of unique variants for the NimbleGen capture panel

# import packages
library(readxl)
library(tidyverse)

# set wd
setwd("~/Desktop/GBM")

# list 3 gbm batches for analysis
gbm_batches <- c("Batch1", "Batch2", "Batch3")

# create unique list of variants for NimbleGen / manual review
total_df <- data.frame()
for ( batch in gbm_batches ) {
    variant_path <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_annotated_variants.xlsx", sep="")
    samples <- excel_sheets(path = variant_path)
    # define tumors
    if (batch == "Batch1" | batch == "Batch2") {
        tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[2])))
    } else {
        tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[1])))
    }
    for ( tumor in tumors ) {
        tumor_df <- c()
        tumor_samples <- samples[grep(tumor, samples)]
        for (name in tumor_samples) {
            # create standardized shorter sample names
            if ( grepl("Re", name) == TRUE ) {
                new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
                #tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
            } else if ( batch == "Batch3") { 
                new_name <- name
                #tumor <- unlist(strsplit(new_name, '-'))[1]
            } else if ( grepl("Re", name) == FALSE ) {
                new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
                #tumor <- unlist(str_split(name, '-'))[2]
            }
            print(c(tumor, new_name))
            # read in file and add to tumor df
            df <- read_excel(variant_path, sheet = name)
            tumor_df <- rbind(tumor_df, df)
        }
        # unique variants only
        uniq_tumor_df <- tumor_df %>% distinct(CHROM, POS, REF, ALT, .keep_all = T)
        # filter based on normal
        uniq_tumor_df <- uniq_tumor_df %>% filter((NORMAL.DP > 30 & NORMAL.DP < 100 & NORMAL.AF < 0.05) | (NORMAL.DP > 100 & NORMAL.AF < 0.025))
        # add to final
        final_df <- data.frame("Tumor" = tumor, uniq_tumor_df)
        total_df <- rbind(total_df, final_df)
    }
}
total_df <- total_df %>% select(-contains('VAF'))

# write df to tsv
write_tsv(total_df, "~/Desktop/GBM/targeted_resequencing/coverage_check/unique_tumor_variants.tsv")

