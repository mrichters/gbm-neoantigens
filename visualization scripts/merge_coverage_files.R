# 08.28.20
# Usage: Rscript <sample name> <vaf annotated tsv> <final tsv file>

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

setwd("~/Desktop/GBM/targeted_resequencing/coverage_check/")

sample_name <- args[1]
exome <- read_tsv(args[2], col_names=T)
new <- read_tsv(args[3], col_names=T)

rc_names <- names(exome)[grepl("AD|AF|DP", names(exome))]
rc <- rc_names[grepl("Tumor|TUMOR", rc_names)]

# change AD to character
exome[rc[grep("AD", rc)]] <- as.character(exome[rc[grep("AD", rc)]])

exome <- exome %>% unite(Variant, CHROM, POS, REF, ALT, sep=':')
new <- new %>% unite(Variant, CHROM, POS, REF, ALT, sep=':')

new <- new %>% unite(AD.combined, AD.REF, AD.ALT, sep=',')
exome_filter <- filter(exome,exome$Variant %in% new$Variant)

for (row in 1:nrow(exome_filter)) {
    var <- as.character(exome_filter[row, "Variant"])
    ts <- filter(new, Variant == var)
    exome_filter[row, rc[grep("DP", rc)]] <- ts$DP
    exome_filter[row, rc[grep("AD", rc)]] <- ts$AD.combined
    exome_filter[row, rc[grep("AF", rc)]] <- ts$VAF
}

# write merged df to file
write.table(exome_filter, paste0(sample_name, ".final.annotated.tsv"), sep='\t', quote=F, row.names=F)
