# 11.4.19
# creating a df of unique variants for the NimbleGen capture panel

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

# create unique list of variants for NimbleGen / manual review
total_shared_value_df <- data.frame()
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
      # add to tumor_df
      df.cbind <- cbind(df[, c("CHROM", "POS", "REF", "ALT", "SYMBOL", "Protein_position", "Amino_acids")], Sample = new_name, Tumor = tumor, Value = 1)
      df.final <- df.cbind %>%
        unite(Variant, CHROM, POS, REF, ALT, sep=':')
      tumor_df <- rbind(tumor_df, df.final)
    }
    # add tumor_df to final total_shared_value_df
    total_shared_value_df <- rbind(total_shared_value_df, tumor_df)
  }
}
# change Sample to character
total_shared_value_df$Sample <- as.character(total_shared_value_df$Sample)
# count unique variants
counts_df <- count(total_shared_value_df[ ,1])
counts_df$x <- as.character(counts_df$x)
# create final df with counts and associated samples, sample[is.na(sample)]
final_counts_df <- data.frame()
for ( row in 1:nrow(counts_df) ) {
  variant <- counts_df[row, 1]
  protein_pos <- filter(total_shared_value_df, Variant == variant)$Protein_position[1]
  amino_acids <- unlist(strsplit(filter(total_shared_value_df, Variant == variant)$Amino_acids[1], '/'))
  if ( is.na(protein_pos) == TRUE || is.na(amino_acids) == TRUE ) {
    protein_change <- '-'
  } else if (length(amino_acids) == 1) {
    protein_change <- paste(amino_acids[1], protein_pos, sep = "")
  } else {
    protein_change <- paste(amino_acids[1], protein_pos, amino_acids[2], sep = "")
  }
  variant_samples <- filter(total_shared_value_df, Variant == variant)$Sample
  genes <- filter(total_shared_value_df, Variant == variant)$SYMBOL[1]
  samples_format <- paste(variant_samples, collapse=", ")
  final_counts_df <- rbind(final_counts_df, cbind(variant, genes, protein_change, counts_df[row, 2], samples_format))
}
# format colnames for final_counts_df
colnames(final_counts_df) <- c("Variant", "Gene", "Protein Position", "Frequency", "Samples")
# save df with separated variant info
separate_df <- separate(final_counts_df, Variant, c("Chrom", "Pos", "Ref", "Alt"))
# write df to excel file
write.xlsx(separate_df, "gbm_unique_variants.xlsx")

################# DO THIS BUT FIND UNIQUE VARIANTS PER TUMOR ###################

# read in cancer type per BrMET
brmet_types <- read_excel("~/Desktop/GBM/gbm_google_drive/brmet_cancer_types.xlsx")

final_tumor_counts_df <- data.frame()
for ( tumor in unique(total_shared_value_df$Tumor) ) {
    if ( grepl("BrMET", tumor) ) { 
        cancer_type <- filter(brmet_types, Tumor == tumor)$`Cancer Type`
        tumor_type <- "BrMET"
    } else if ( grepl("Re", tumor) ) {
        tumor_type <- "Recurrent GBM"
        cancer_type <- "Recurrent GBM"
    } else {
        tumor_type <- "GBM"
        cancer_type <- "GBM"
    }
    # count unique variants
    tumor_total_df <- filter(total_shared_value_df, Tumor == tumor)
    counts_df <- count(tumor_total_df[ ,1])
    colnames(counts_df) <= c("Variant", "Count")
    final_tumor_counts_df <- rbind(final_tumor_counts_df, data.frame("Tumor" = tumor, "Tumor Type" = tumor_type, "Cancer Type" = cancer_type, "Unique_Counts" = nrow(counts_df)))
}
# order by cancer type
final_tumor_counts_df$Cancer.Type <- factor(final_tumor_counts_df$Cancer.Type, levels = c("Breast Cancer", "Melanoma", "NSCLC", "SCLC", "GBM", "Recurrent GBM"))
final_tumor_counts_df <- final_tumor_counts_df[order(final_tumor_counts_df$Cancer.Type), ]

colnames(final_tumor_counts_df) <- gsub('\\.', " ", colnames(final_tumor_counts_df))
# plot
ggplot(final_tumor_counts_df, aes(x=`Tumor Type`, y=Unique_Counts)) + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + geom_boxplot(notch=FALSE, outlier.shape = NA,fill="grey90") + geom_jitter(position=position_jitter(0.2),aes(color=`Cancer Type`), size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=12),axis.title.x = element_blank()) + ylab("Variant Count") + scale_color_manual(values = pal_jama("default", alpha=0.9)(7)[1:7]) + coord_trans(y="log2") + scale_y_continuous(breaks=c(100,500,1000,2000,4000,6000))


################# FORMATTING LIFTOVER VARIANTS from final file ^^ ####################
setwd("~/Desktop/GBM/targeted_resequencing/")

previous_output <- as.data.frame(read_excel("gbm_unique_variants_b38.xlsx"))
previous_output.u <- unite(previous_output, "Coordinate", Chrom:Pos, sep = ":")
duplicates <- data.frame("Duplicated_coordinates" = sort(previous_output.u$Coordinate[duplicated(previous_output.u$Coordinate)]))
unique_all_b38 <- data.frame("Coordinates" = unique(previous_output.u$Coordinate))

failed_liftover <- as.data.frame(read_tsv("liftover_fails.vcf", col_names=F))
failed_liftover.u <- unite(failed_liftover, "Coordinate", X1:X2, sep = ":")

failed_vector <- as.vector(failed_liftover.u$Coordinate)
failed_vector_unique <- unique(failed_vector) 

failed_present <- c()
for ( item in unique_all_b38$Coordinates ) {
  if ( item %in% failed_vector ) {
    failed_present <- c(failed_present, item)
  }
}

failed_absent <- c()
for ( item in failed_vector_unique ) {
  if ( is.element(item,failed_present) == F )
    failed_absent <- c(failed_absent, item)
} 
  
bed_format <- read_tsv("Megan/unique_variants_b38_coordinates_vcf_to_bed.bed", col_names = F)
vcf_format <- read_tsv("Megan/unique_variants_b38_coordinates_bed_to_vcf.vcf", col_names = F)
vcf_format.u <- unite(vcf_format, "Coordinate", X1:X2, sep = ":")
vcf_format.u_vector <- as.vector(vcf_format.u$Coordinate)

for ( item in unique_all_b38$Coordinates ) {
  if ( is.element(item, vcf_format.u_vector) == F ) {
    print(item)
  }
}

final_df <- data.frame()
for ( row in 1:nrow(previous_output.u) ) {
  coor <- previous_output.u[row, "Coordinate"]
  if ( coor %in% failed_vector_unique ) {
    final_df_line <- data.frame(previous_output.u[row, ], "FILTER" = "Fail")
  } else {
    final_df_line <- data.frame(previous_output.u[row, ], "FILTER" = "Pass")
  }
    final_df <- rbind(final_df, final_df_line)
}
final_df_sep <- separate(final_df, Coordinate, c("Chrom", "Pos")) 
   
write.xlsx(final_df_sep, "gbm_unique_variants_b38_filter.xlsx")

