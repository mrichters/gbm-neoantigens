#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# load packages and set wd
library(dplyr)
library(tidyr)

#set wd
setwd(args[1])  

# Read in table made with VariantsToTable command
normal_calls <- read.delim(args[2], header = TRUE, col.names = c("CHROM", "POS", "GT", "AD", "DP"))

# Filter for positions with at least 20x coverage
normal_calls <- subset.data.frame(normal_calls, DP >= 20)

# Split AD into REF and ALT 
AD_REF <- sapply(normal_calls$AD, function(x) return(unlist(strsplit(x, ","))[1]))
AD_ALT <- sapply(normal_calls$AD, function(x) return(unlist(strsplit(x, ","))[2]))
normal_calls$REF <- as.numeric(AD_REF)
normal_calls$ALT <- as.numeric(AD_ALT)

# Calculate VAFs, allele depth/total depth
normal_calls$VAF <- normal_calls$ALT/normal_calls$DP

# Filter for heterozygous posistions, i.e. positions with VAF between .4 and .6
normal_calls <- subset.data.frame(normal_calls, VAF >= .4); normal_calls <- subset.data.frame(normal_calls, VAF <= .6)

# Create table containing info for heterozygous SNPs, will use later to plot results
write.table(normal_calls, file="het.snps.table", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Get list of heterozygous positions to use with bam-readcount
heterozygous.positions <- unique(normal_calls[c("CHROM", "POS")])
write.table(heterozygous.positions[ ,c("CHROM", "POS", "POS")], file = "het.positions.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Exit R, you do not need to save workspace 
q()  
