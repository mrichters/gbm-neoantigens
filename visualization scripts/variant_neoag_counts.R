# 12.13.19
# counts

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
setwd("~/Desktop/GBM")

# list 3 gbm batches for analysis
gbm_batches <- c("Batch1", "Batch2", "Batch3")

# counts
counts_df <- data.frame()
sample_names <- c() #shorter sample names
for ( batch in gbm_batches ) {
  variant_path <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_annotated_variants.xlsx", sep="")
  samples <- excel_sheets(path = variant_path)
  tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[2])))
  for ( tumor in tumors ) {
    tumor_df <- c()
    for ( name in samples[grep(tumor, samples)] ) {
      # read in correct batch variant excel file
      df <- as.data.frame(read_excel(variant_path, sheet = name))
      # create standardized shorter sample names
      if ( grepl("Re", name) == TRUE ) {
        new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
        tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
      } else if ( batch == "Batch3") { 
          new_name <- mgsub(name, c("19-", "Tumor", "_"), c("", "", "-"))
          tumor <- unlist(strsplit(new_name, '-'))[1]
      } else if ( grepl("Re", name) == FALSE ) {
        new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
        tumor <- unlist(str_split(name, '-'))[2]
      } 
      # add to sample_names df
      sample_names <- c(sample_names, new_name)
      df.final <- data.frame("Tumor" = tumor, "Sample" = new_name, "Variant Count" = nrow(df))
      tumor_df <- rbind(tumor_df, df.final)
    }
    # add tumor_df to final total_shared_value_df
    counts_df <- rbind(counts_df, tumor_df)
  }
}
# deal with factors to get samples in the correct order
counts_df$Sample <- factor(counts_df$Sample, levels = sort(as.character(counts_df$Sample)))
counts_df$Tumor <- factor(counts_df$Tumor, levels = sort(unique(as.character(counts_df$Tumor))))
# order new df
counts_df.ordered <- counts_df[order(counts_df$Sample), ]
# customize colors with ggsci
my_pal <- c(pal_d3("category10", alpha = 1)(10), pal_aaas("default", alpha = 1)(7), pal_d3("category20c", alpha = 1)(10))
counts_df_065 <- filter(counts_df, Tumor == "Re.GBM065")
counts_df_065$Sample <- c("Re.GBM056-primary", "Re.GBM065-1","Re.GBM065-2", "Re.GBM065-3", "Re.GBM065-4", "Re.GBM065-5")
counts_df_065$`Tumor Type` <- c("Primary", "Recurrent", "Recurrent", "Recurrent", "Recurrent", "Recurrent")
colnames(counts_df_065) <- gsub('\\.', " ", colnames(counts_df_065))
# bar plot of counts
# BRMETs + hypermutator
ggplot(counts_df_065, aes(x=Sample, y=`Variant Count`, fill = `Tumor Type`)) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + scale_fill_manual(values = c("darkgreen", "grey40"))
# GBMs
ggplot(filter(counts_df[grepl("GBM", counts_df$Sample) & !grepl("Re.GBM065", counts_df$Sample), ]), aes(x=Sample, y=Variant.Count, fill = Tumor)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), ) + scale_fill_manual(values = my_pal) + ggtitle("Primary and Recurrent GBM Variant Counts")

# 1.15.20
# only hypermutator variant counts - primary and recurrent
counts_df <- data.frame()
variant_path <- "~/Desktop/GBM/GBM065/GBM065_shared.xlsx"
samples <- excel_sheets(path = variant_path)
tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[1])))
for ( tumor in tumors ) {
  tumor_df <- c()
  for ( name in samples[grep(tumor, samples)] ) {
    # read in correct batch variant excel file
    df <- as.data.frame(read_excel(variant_path, sheet = name))
    df.final <- data.frame("Tumor" = tumor, "Sample" = name, "Variant Count" = nrow(df))
    tumor_df <- rbind(tumor_df, df.final)
  }
  counts_df <- rbind(counts_df, tumor_df)
}

counts_df$Tumor_type <- c("Primary", "Recurrent", "Recurrent", "Recurrent")
counts_df$Sample <- c("Primary", "Recurrence-1", "Recurrence-2", "Recurrence-3")

ggplot(counts_df, aes(x=Sample, y=Variant.Count, fill = Tumor_type)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position = "none", panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), text = element_text(size=15)) + ylab("Variant Count") + scale_fill_manual(values = c("darkolivegreen", "grey30"))

neoantigen_df <- as.data.frame(read_excel("~/Desktop/GBM/neoantigen_forming_variants_proportion.xlsx"))
neoantigen_df_065 <- neoantigen_df[grep("GBM065", neoantigen_df$Sample), ]
neoantigen_df_065$Tumor_type <- "Recurrent"
neoantigen_df_065$Sample <- c("Recurrence-1", "Recurrence-2", "Recurrence-3")

ggplot(neoantigen_df_065, aes(x=Sample, y=Total_Neoantigens)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position = "none", panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), text = element_text(size=15)) + scale_color_manual(values = "grey50") + ylab("Neoantigen Count") + xlab("Sample")

# 1.15.20
# creating boxplots of brmet variants counts vs gbm variant counts
tumor_type_vector <- c()
for ( row in 1:nrow(counts_df) ) {
  tumor <- counts_df[row, "Tumor"]
  if ( grepl("BrMET", tumor) ) {
    tumor_type <- "BrMET"
  } else if ( grepl("Re", tumor) ) {
    tumory_type <- "Recurrent GBM"
  } else {
    tumor_type <- "GBM"
  }
  tumor_type_vector <- c(tumor_type_vector, tumor_type)
}
counts_df$Tumor_Type <- tumor_type_vector

ggplot(counts_df, aes(x=Tumor_Type, y=Variant.Count)) + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + geom_boxplot(notch=FALSE, outlier.shape = NA, width=0.3) + geom_jitter(position=position_jitter(0.15),aes(color=Tumor_Type), size=2.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=15), legend.position = "none") + ylab("Variant Count") + xlab("Tumor Type")+ scale_color_manual(values = c("indianred4", "grey30")) 
#[!grepl("Re.GBM065", counts_df$Sample), ]
