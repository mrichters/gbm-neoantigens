#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Set working directory
setwd(args[1])

# Read bam-readcount output into R
tumor_VAFs <- read.delim(args[2], header = FALSE, fill = TRUE, col.names = c("CHROM", "POS", "REF", "TUMOR_DP", "5", "count_A", "count_C", "count_G", "count_T", "11", "12", "13", "14"))

# Select columns relevant to our calculation
tumor_VAFs <- select(tumor_VAFs, CHROM, POS, TUMOR_DP, count_A, count_C, count_G, count_T)

# Filter down to positions with at least 20x coverage
tumor_VAFs <- subset.data.frame(tumor_VAFs, TUMOR_DP >= 20)

# Spilt columns with allele info into separate columns with allele depth
tumor_VAFs <- tumor_VAFs %>% separate(count_A, c("A", "A_AD"), extra = "drop"); tumor_VAFs <- tumor_VAFs %>% separate(count_C, c("C", "C_AD"), extra = "drop"); tumor_VAFs <- tumor_VAFs %>% separate(count_G, c("G", "G_AD"), extra = "drop"); tumor_VAFs <- tumor_VAFs %>% separate(count_T, c("T", "T_AD"), extra = "drop")

# Again, select columns relevant to VAF calculation
tumor_VAFs <- select(tumor_VAFs, CHROM, POS, TUMOR_DP, A_AD, C_AD, G_AD, T_AD)
 
# Make sure allele depth columns are numeric
tumor_VAFs$A_AD <- as.numeric(tumor_VAFs$A_AD); tumor_VAFs$C_AD <- as.numeric(tumor_VAFs$C_AD); tumor_VAFs$G_AD <- as.numeric(tumor_VAFs$G_AD); tumor_VAFs$T_AD <- as.numeric(tumor_VAFs$T_AD)

# Calculate VAFs
tumor_VAFs$VAF_A <- tumor_VAFs$A_AD/tumor_VAFs$TUMOR_DP; tumor_VAFs$VAF_C <- tumor_VAFs$C_AD/tumor_VAFs$TUMOR_DP; tumor_VAFs$VAF_G <- tumor_VAFs$G_AD/tumor_VAFs$TUMOR_DP; tumor_VAFs$VAF_T <- tumor_VAFs$T_AD/tumor_VAFs$TUMOR_DP

print(tumor_VAFs)

# Split into separate files for each VAF, remerge so we have one column with all tumor VAFs
VAF_A <- select(tumor_VAFs, CHROM, POS, TUMOR_DP, A_AD, VAF_A); VAF_C <- select(tumor_VAFs, CHROM, POS, TUMOR_DP, C_AD, VAF_C); VAF_G <- select(tumor_VAFs, CHROM, POS, TUMOR_DP, G_AD, VAF_G); VAF_T <- select(tumor_VAFs, CHROM, POS, TUMOR_DP, T_AD, VAF_T)
colnames(VAF_A) <- c("CHROM", "POS", "TUMOR_DP", "AD", "VAF"); colnames(VAF_C) <- c("CHROM", "POS", "TUMOR_DP", "AD", "VAF"); colnames(VAF_G) <- c("CHROM", "POS", "TUMOR_DP", "AD", "VAF"); colnames(VAF_T)<- c("CHROM", "POS", "TUMOR_DP", "AD", "VAF")
tumor_VAFs <- merge(VAF_A, VAF_C, all = TRUE); tumor_VAFs <- merge(tumor_VAFs, VAF_G, all = TRUE); tumor_VAFs <- merge(tumor_VAFs, VAF_T, all = TRUE)

# Filter out VAFs where no reads mapped (i.e VAFs = 0)
tumor_VAFs <- subset.data.frame(tumor_VAFs, VAF > 0)

# Create table for tumor VAFs
write.table(tumor_VAFs, file=args[3], sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Read in previously made table with germline VAFs
#germline_VAFs <- read.delim(args[4], header = TRUE, col.names = c("CHROM", "POS", "GT", "AD", "DP", "VAF"))

# Plot VAFs
#pdf(args[5], width = 11.5, height = 2.75)
#VAF_plot <- ggplot() + geom_point(data = germline_VAFs, aes(POS,VAF), color="blue", size = .75) + geom_point(data = tumor_VAFs, aes(POS,VAF), color="green", size = .75) + ylim(0,1) + xlab(args[6]) + ylab("VAF") + ggtitle(args[7])
#plot(VAF_plot)
#dev.off()
