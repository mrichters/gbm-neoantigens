# GGpairs VAF scatter plots #
# input is excel file with annotated variants #

library(GGally)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggsci)

# functions
filter_variants <- function(variant_df) {
    filtered_df <- data.frame()
    for (row in 1:nrow(variant_df)) {
        caller_set <- unlist(strsplit(variant_df[row, "set"], "-"))
        if (length(caller_set) > 1) {
            filtered_df <- rbind(filtered_df, variant_df[row, ])
        }
    }
    return(filtered_df)
}

setwd("~/Desktop/GBM\ Figures/")

variant_path <- "~/Desktop/GBM/gbm_google_drive/merged_annotated_variants.xlsx"
samples <- excel_sheets(path = variant_path)
tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[1])))

for ( tumor in "GBM065" ) {
    tumor_df <- c()
    sample_names <- c()
    for ( name in samples[grep(tumor, samples)] ) {
        # read in correct batch variant excel file
        df <- as.data.frame(read_excel(variant_path, sheet = name))
        # if want to filter results to only include those called by 2 variant callers
        #df <- filter_variants(df)
        new_name <- name
        df.subset <- df %>% 
            select(CHROM, POS, REF, ALT, `TUMOR.AF`)
        df.subset$`TUMOR.AF` <- as.numeric(df.subset$`TUMOR.AF`)
        names(df.subset)[5] <- "VAF"
        df.subset.cbind <- cbind(df.subset, Sample = new_name)
        df.final <- df.subset.cbind %>%
            unite(Variant, CHROM, POS, REF, ALT, sep=':')
        tumor_df <- rbind(tumor_df, df.final)
    }
    df.spread <- tumor_df %>% 
        spread(Sample, VAF)
    df.spread[is.na(df.spread)] <- 0
    df.ggpairs <- df.spread[,-c(1)]
    # add dots color
    df.ggpairs$Gene <- "Other"
    driver_coordinates <- c("chr17:7673782:T:C", "chr17:7675994:C:T", "chr2:208248388:C:T", "chr2:47806213:C:T")
    driver_indices <- sapply(driver_coordinates, function(x) return(match(x,df.spread$Variant)))
    df.ggpairs[driver_indices,]$Gene <- c("TP53 R280G", "TP53 T125T", "IDH1 R132H", "MSH6 T1219I")
    # factor and order
    df.ggpairs$Gene <- factor(df.ggpairs$Gene, levels = c("Other", "TP53 R280G", "TP53 T125T", "IDH1 R132H", "MSH6 T1219I"))
    df.ggpairs <- df.ggpairs[order(df.ggpairs$Gene),]
    my_palette <- c("grey30", pal_locuszoom()(5)[c(1:4)])
    pm <- ggpairs(df.ggpairs, diag=list(continuous="barDiag"), axisLabels='show', upper='blank', columns=1:5, mapping = ggplot2::aes(color=Gene), legend=c(2,1))
    for(i in 2:pm$nrow) {
        for(j in 1:(i-1)) {
            pm[i,j] <- pm[i,j] +
                scale_x_continuous(limits = c(0, 1)) +
                scale_y_continuous(limits = c(0, 1)) +
                scale_color_manual(values = my_palette) +
                theme(legend.title = element_blank())
        }
    }
    for(i in 1:5) {
        for(j in 1:5) {
            pm[i,j] <- pm[i,j] +
                scale_fill_manual(values = my_palette)
        }
    }
    # scatter plot object #
    pm
    #pdf(paste(tumor, "_vaf_scatter_discrete_unfiltered_tr.pdf", sep = ""), width = 10, height = 8)
    #print(pm)
    #dev.off()
}
