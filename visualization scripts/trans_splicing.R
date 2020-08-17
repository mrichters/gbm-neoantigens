# 8.29.19

# import packages
library(readxl)
library(tidyverse)
library(ggsci)
library(plyr); library(dplyr)
library(gridExtra)
library(gtable)
library(scales)
library(GenVisR)

#set wd
setwd("~/Desktop/GBM_plots")

splicing_gene_set <- c("AGGF1", "C9orf78", "CCAR1", "CD2BP2", "CDC5L", "CDK11A", "CDK12", "CELF4", "CFAP20", "CLK4", "CWC22", "DDX17", "DDX18", "DDX20", "DDX23", "DDX26B", "DDX27", "DDX3X", "DDX41", "DDX5", "DDX50", "DHX16", "DHX35", "DHX36", "DHX9", "EEF1A1", "EFTUD2", "EIF2S2", "ELAVL1", "ELAVL2", "ELAVL4", "FAM58A", "FRA10AC1", "FUBP1", "FUBP3", "GPATCH8", "HNRNPCL1", "HNRNPD", "HNRNPDL", "HNRNPH3", "HNRNPK", "HNRNPL", "IGF2BP3", "INTS4", "INTS7", "KIAA1429", "KIN", "MBNL2", "MOV10", "NCBP1", "NELFE", "NOVA1", "NSRP1", "PABPC1", "PCBP1", "PCBP2", "PCBP3", "PHF5A", "PLRG1", "PPIG", "PPIL1", "PPIL4", "PRPF3", "PRPF38B", "PRPF39", "PRPF40B", "PRPF4B", "PSIP1", "QKI", "RALYL", "RBBP6", "RBM10", "RBM15B", "RBM25", "RBM26", "RBM27", "RBM7", "RBM8A", "RBMX", "RBMX2", "RNF20", "SF1", "SF3B1", "SF3B2", "SF3B3", "SKIV2L2", "SNRNP200", "SNRNP35", "SNRNP48", "SNRPD3", "SNRPN", "SPEN", "SRSF2", "SRSF5", "SYNCRIP", "TCERG1", "THOC5", "THOC6", "THOC7", "THRAP3", "TIA1", "TIAL1", "TNPO1", "TRIM24", "TTC14", "U2AF1", "U2AF2", "U2SURP", "WBP11", "WBP4", "ZC3H13", "ZC3H18", "ZC3H4", "ZCCHC8", "ZCRB1", "ZMYM3", "ZNF131", "ZNF207", "ZRSR2")

# list all sheets in excel file
twck_samples <- excel_sheets(path="~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx")
twck_tumors <- unique(sapply(twck_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
htc_samples <- excel_sheets(path="~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx")
htc_tumors <- unique(sapply(htc_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
ribas_samples <- excel_sheets(path="~/Desktop/GBM/ribas_new_annotated_tsvs/Ribas_annotated_variants.xlsx")

#htc
variants <- "~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx"
#twck
variants <- "~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx"
# ribas melanoma
variants <- "~/Desktop/GBM/ribas_new_annotated_tsvs/Ribas_annotated_variants.xlsx"

# run twice - once with HTC and once with TWCK
total_df <- c()
mut_df <- c()
samples_df <- c()
for ( tumor in htc_tumors ) {
  tumor_df <- c()
  for ( name in htc_samples[grep(tumor, htc_samples)] ) {
    gbm_data <- as.data.frame(read_excel(variants, sheet = name))
    if ( grepl("Re", name) == FALSE ) {
      new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
      tumor <- unlist(str_split(name, '-'))[2]
      if ( grepl("BrMET", name) ) {
        tumor_type <- "Brain Metastasis"
      } else {
        tumor_type <- "Primary GBM"
      }
    } else { 
      new_name <- paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
      tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
      tumor_type <- "Recurrent GBM"
    } 
    gbm_mutation_data <- gbm_data[,c("SYMBOL", "Consequence", "HGVSp")]
    gbm_mutation_data$Consequence <- unlist(lapply(strsplit(gbm_mutation_data$Consequence, "&"), `[[`, 1))
    gbm_mutation_data[is.na(gbm_mutation_data)] <- 0
    gbm_mutation_data_filter <- filter(gbm_mutation_data, HGVSp != 0 & SYMBOL %in% splicing_gene_set & Consequence != "synonymous_variant")
    gbm_mutation_data_filter$HGVSp <- unlist(lapply(strsplit(gbm_mutation_data_filter$HGVSp, ":"), `[[`, 2))
    if ( nrow(gbm_mutation_data_filter) > 0 ) {
      gbm_mutation_data_filter <- cbind("sample" = new_name, gbm_mutation_data_filter)
      colnames(gbm_mutation_data_filter) <- c("sample", "gene", "variant_class", "amino.acid.change")
      tumor_df <- rbind(tumor_df, gbm_mutation_data_filter)
    }
    mut_df_line <- cbind(as.data.frame(new_name), as.data.frame(nrow(gbm_mutation_data)))
    colnames(mut_df_line) <- c("sample", "mut_burden")
    mut_df <- rbind(mut_df, mut_df_line)
    samples_line <- data.frame(new_name, "Tumor Type", tumor_type)
    colnames(samples_line) <- c("sample", "variable", "value")
    samples_df <- rbind(samples_df, samples_line)
  }
  total_df <- rbind(total_df, tumor_df)
}
# create sorted sample order
sample_order <- sort(as.vector(mut_df$sample))
# change sample order
total_df$sample <- as.character(total_df$sample)
total_df.sort <- total_df[order(total_df$sample),]
mut_df$sample <- as.character(mut_df$sample)
mut_df.sort <- mut_df[order(mut_df$sample), ]
# most frequently mutated genes
table_genes <- as.data.frame(table(total_df$gene))
table.sort <- table_genes[order(-table_genes$Freq), ]
top_genes_list <- as.vector(table.sort$Var1[1:50])
total_df_genes <- filter(total_df, gene %in% top_genes_list & variant_class != "protein_altering_variant")
mut_df.subset <- filter(mut_df, sample %in% total_df_genes$sample)
mut_df.subset.sort <- mut_df.subset[order(mut_df.subset$sample), ]
samples_df.subset <- filter(samples_df, sample %in% total_df_genes$sample)
samples_df.subset.sort <- samples_df.subset[order(samples_df.subset$sample), ]
sample_ordering_waterfall <- sort(as.vector(mut_df.subset$sample))
# Create a vector to save mutation priority order for plotting
#mutation_priority <- as.character(unique(total_df_genes$variant_class))
mutation_priority <- c("frameshift_variant", "missense_variant", "stop_gained")
mutationColours <- c("frameshift_variant"='#A80100', "stop_gained"='#4f00A8', "missense_variant"='#009933')
# Create an trans splicing waterfall plot - gbm
waterfall(total_df, fileType = "Custom", variant_class_order = mutation_priority, mainXlabel = TRUE, mutBurden = mut_df.subset.sort, mainPalette=mutationColours, clinData = samples_df.subset.sort, clinLegCol=1, clinVarCol=c('Brain Metastasis'='#c2ed67', 'Primary GBM'='#E63A27', 'Recurrent GBM'='#e69127'), clinVarOrder=c('Brain Metastasis', 'Primary GBM', 'Recurrent GBM'), section_heights=c(1, 6, 0.5), sampOrder = sample_ordering_waterfall)


#run for ribas melanoma
total_df <- c()
mut_df <- c()
for ( name in ribas_samples ) {
  new_name <- name
  gbm_data <- as.data.frame(read_excel(variants, sheet = name))
  gbm_mutation_data <- gbm_data[,c("SYMBOL", "Consequence", "HGVSp")]
  gbm_mutation_data$Consequence <- unlist(lapply(strsplit(gbm_mutation_data$Consequence, "&"), `[[`, 1))
  gbm_mutation_data[is.na(gbm_mutation_data)] <- 0
  gbm_mutation_data_filter <- filter(gbm_mutation_data, HGVSp != 0 & SYMBOL %in% splicing_gene_set & Consequence != "synonymous_variant")
  gbm_mutation_data_filter$HGVSp <- unlist(lapply(strsplit(gbm_mutation_data_filter$HGVSp, ":"), `[[`, 2))
  if ( nrow(gbm_mutation_data_filter) > 0 ) {
    gbm_mutation_data_filter <- cbind("sample" = new_name, gbm_mutation_data_filter)
    colnames(gbm_mutation_data_filter) <- c("sample", "gene", "variant_class", "amino.acid.change")
    total_df <- rbind(total_df, gbm_mutation_data_filter)
  }
  mut_df_line <- cbind(as.data.frame(new_name), as.data.frame(nrow(gbm_mutation_data)))
  colnames(mut_df_line) <- c("sample", "mut_burden")
  mut_df <- rbind(mut_df, mut_df_line)
}
# create sorted sample order
sample_order <- sort(as.vector(mut_df$sample))
# change sample order
total_df$sample <- as.character(total_df$sample)
total_df.sort <- total_df[order(total_df$sample),]
mut_df$sample <- as.character(mut_df$sample)
mut_df.sort <- mut_df[order(mut_df$sample), ]
mut_df.subset <- filter(mut_df, sample %in% total_df$sample)
mut_df.subset.sort <- mut_df.subset[order(mut_df.subset$sample), ]
sample_ordering_waterfall <- sort(as.vector(mut_df.subset$sample))
# Create a vector to save mutation priority order for plotting
#mutation_priority <- as.character(unique(total_df_genes$variant_class))
mutation_priority <- c("frameshift_variant", "missense_variant", "stop_gained")
mutationColours <- c("frameshift_variant"='#A80100', "stop_gained"='#4f00A8', "missense_variant"='#009933')

# Create an trans splicing waterfall plot - melanoma
waterfall(total_df, fileType = "Custom", variant_class_order = mutation_priority, mainXlabel = TRUE, mutBurden = mut_df.subset.sort, mainPalette=mutationColours, sampOrder = sample_ordering_waterfall)
