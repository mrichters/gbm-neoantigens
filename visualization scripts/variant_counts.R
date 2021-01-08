# 11.4.19
# variant counts

# import packages
library(readxl)
library(tidyverse)
library(ggsci)
library(plyr); library(dplyr)
library(GGally)
library(UpSetR)
library(gridExtra)
library(gtable)
library(scales)
library(openxlsx)
library(mgsub)

# set wd
setwd("~/Desktop/GBM Figures/")

# samples / tumors in cohort
variant_path <- paste0("~/Desktop/GBM/gbm_google_drive/", "merged_annotated_variants.xlsx")
samples <- excel_sheets(path = variant_path)
tumors <- unique(sapply(samples, function(x) return(unlist(strsplit(x, '-'))[1])))

# create unique list of variants per tumor for variant counts boxplot
total_tumor_variants <- data.frame()
# loop over tumors
for ( tumor in tumors ) {
    tumor_df <- c()
    # loop over samples in tumor
    for ( name in samples[grep(tumor, samples)] ) {
      # read in correct batch variant excel file
      df <- as.data.frame(read_excel(variant_path, sheet = name))
      # add sample name, tumor, and count info
      df.cbind <- cbind(df[, c("CHROM", "POS", "REF", "ALT", "SYMBOL", "Protein_position", "Amino_acids")], Sample = name, Tumor = tumor, Value = 1)
      # unite variant info for comparison to other variants in tumor
      df.final <- df.cbind %>% unite(Variant, CHROM, POS, REF, ALT, sep=':')
      # add to tumor df
      tumor_df <- rbind(tumor_df, df.final)
    }
    # add tumor_df to final total_shared_value_df
    total_tumor_variants <- rbind(total_tumor_variants, tumor_df)
}

## FOR VARIANT COUNTS ## 
# read in cancer type per BrMET
brmet_types <- read_excel("~/Desktop/GBM/gbm_google_drive/brmet_cancer_types.xlsx")

# assign tumor / cancer types
final_tumor_counts_df <- data.frame()
for ( tumor in unique(total_tumor_variants$Tumor) ) {
    if ( grepl("BrMET", tumor) ) { 
        cancer_type <- filter(brmet_types, Tumor == tumor)$`Cancer Type`
        tumor_type <- "BrMET"
    } else if ( grepl("Re", tumor) ) {
        tumor_type <- "GBM"
        cancer_type <- "Recurrent GBM"
    } else {
        tumor_type <- "GBM"
        cancer_type <- "GBM"
    }
    # count unique variants
    tumor_total_df <- filter(total_tumor_variants, Tumor == tumor)
    counts_df <- count(tumor_total_df[ ,1])
    colnames(counts_df) <= c("Variant", "Count")
    final_tumor_counts_df <- rbind(final_tumor_counts_df, data.frame("Tumor" = tumor, "Tumor Type" = tumor_type, "Cancer Type" = cancer_type, "Unique_Counts" = nrow(counts_df)))
}

# order by cancer type
final_tumor_counts_df$Cancer.Type <- factor(final_tumor_counts_df$Cancer.Type, levels = c("Breast Cancer", "Melanoma", "NSCLC", "SCLC", "GBM", "Recurrent GBM"))
final_tumor_counts_df <- final_tumor_counts_df[order(final_tumor_counts_df$Cancer.Type), ]
colnames(final_tumor_counts_df) <- gsub('\\.', " ", colnames(final_tumor_counts_df))

# color palette
mypal <- c("#E64B35FF", "#7E6148FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF")

# plot
variant_counts_plot <- ggplot(final_tumor_counts_df, aes(x=`Tumor Type`, y=Unique_Counts)) + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + geom_violin(fill = "grey90") + geom_jitter(position=position_jitter(0.1),aes(color=`Cancer Type`), size=2) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8),axis.title.x = element_blank(), legend.title = element_blank()) + ylab("Variant Count") + scale_color_manual(values = c("#E64B35FF", "#7E6148FF", "#00A087FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF")) + coord_trans(y="log2") + scale_y_continuous(breaks=c(100,500,1000,2000,4000,6000))

# geom_boxplot(notch=FALSE, outlier.shape = NA,fill="grey90")

