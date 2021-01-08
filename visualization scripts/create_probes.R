# 11.4.19
# creating a df of unique variants for the NimbleGen capture panel
# and post - targeted seq

# import packages
library(readxl)
library(tidyverse)
library(plyr); library(dplyr)
library(xlsx)
library(arsenal)

# set wd
setwd("~/Desktop/GBM Figures/")

# samples / tumors in cohort
variant_path <- paste0("~/Desktop/GBM/gbm_google_drive/merged_annotated_variants.xlsx")
samples <- excel_sheets(path = variant_path)
tumors <- unique(sapply(samples, function(x) return(unlist(strsplit(x, '-'))[1])))

# create unique list of variants per tumor for NimbleGen / manual review 
total_shared_value_df <- data.frame()
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
    total_shared_value_df <- rbind(total_shared_value_df, tumor_df)
}

# count unique variants
counts_df <- count(total_shared_value_df[ ,1])
counts_df$x <- as.character(counts_df$x)
# create final df with counts and associated samples, sample[is.na(sample)]
final_counts_df <- data.frame()
# loop over lines in counts_df
for ( row in 1:nrow(counts_df) ) {
    # define variables
    variant <- counts_df[row, 1]
    protein_pos <- filter(total_shared_value_df, Variant == variant)$Protein_position[1]
    amino_acids <- unlist(strsplit(filter(total_shared_value_df, Variant == variant)$Amino_acids[1], '/'))
    # create AA change column
    if ( is.na(protein_pos) == TRUE || is.na(amino_acids) == TRUE ) {
        protein_change <- '-'
    } else if (length(amino_acids) == 1) {
        protein_change <- paste(amino_acids[1], protein_pos, sep = "")
    } else {
        protein_change <- paste(amino_acids[1], protein_pos, amino_acids[2], sep = "")
    }
    # get sample that contains variant 
    variant_samples <- filter(total_shared_value_df, Variant == variant)$Sample
    samples_format <- paste(variant_samples, collapse=", ")
    # gene that variant is in
    genes <- filter(total_shared_value_df, Variant == variant)$SYMBOL[1]
    # add row to final_counts_df - 1 line per unique variant
    final_counts_df <- rbind(final_counts_df, cbind(variant, genes, protein_change, counts_df[row, 2], samples_format))
}
# format colnames for final_counts_df
colnames(final_counts_df) <- c("Variant", "Gene", "Protein Position", "Frequency", "Samples")
# save df with separated variant info
separate_df <- separate(final_counts_df, Variant, c("Chrom", "Pos", "Ref", "Alt"))
# write df to excel file
write.xlsx(separate_df, "ts_validated_gbm_unique_variants.xlsx")

### compare pre- post- targeted seq unique variant lists ###
dir <- "~/Desktop/GBM/gbm_google_drive/targeted_reseq_variants/"
pre <- as.data.frame(read_excel(paste0(dir,"gbm_unique_variants.xlsx")))
post <- as.data.frame(read_excel(paste0(dir, "ts_validated_gbm_unique_variants.xlsx")))
merged <- merge(pre[,c(1:6)], post[,c(1:6)])
# 10,254

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