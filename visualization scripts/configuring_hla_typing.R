# 11.15.19
# configure clinical hla typing results for pvacseq input
# class I done

# import packages
library(tidyverse)
library(readxl)
library(openxlsx)
library(ggsci)

# set wd
setwd("~/Desktop/GBM/hla_typing/batch_final_histogenetics_results/")

### INPUT DATA ###
batch3_typing <- as.data.frame(read_excel("Batch3_Histogenetics_Report.xlsx"), col_names = T)
hla_key <- as.data.frame(read_tsv("../hla_alleles_key.format.tsv"), header = T)


### CLASS I ###
classI <- batch3_typing[, c(1:7)]
alleles <- colnames(classI)[2:ncol(classI)]
final_df <- data.frame()
for ( row in 1:nrow(classI) ) {
  sample <- classI[row, 1]
  df.format <- data.frame()
  for ( item in alleles ) {
    allele_name <- substr(classI[row, item], 1, 5)
    allele_format <- data.frame(sample, paste("HLA-", substr(item, 1, 1), "*", allele_name, sep = ""))
    colnames(allele_format) <- c("Sample", item)
    if ( nrow(df.format) == 0 ) {
      df.format <- allele_format
    } else {
      df.format <- merge(df.format, allele_format)
    }
  }
  final_df <- rbind(final_df, df.format)
}
final_df.format <- unite(final_df, "HLA Class I", A1:C2, remove=T, sep=',')
write_tsv(final_df.format, "../batch3_classI_hlatypes.tsv")


### CLASS II ###
classII <- batch3_typing[, c(1,8:ncol(batch3_typing))]
alleles <- colnames(classII)[2:ncol(classII)]
final_df <- data.frame()
for ( row in 1:nrow(classII) ) {
  sample <- classII[row, 1]
  df.format <- data.frame()
  counter = 0
  for ( item in alleles ) {
    all_alleles <- classII[row, item]
    if ( is.na(all_alleles) == F ) {
      all_alleles.split <- unlist(strsplit(all_alleles, "/"))
      for ( allele in all_alleles.split ) {
        allele_name <- paste(unlist(strsplit(allele, ":"))[1], unlist(strsplit(allele, ":"))[2], sep=":")
        allele_format <- data.frame(sample, paste(substr(item, 1, nchar(item)-1), "*", allele_name, sep = ""))
        colnames(allele_format) <- c("Sample", paste("allele", counter, sep="_"))
        counter <- counter + 1
        if ( nrow(df.format) == 0 ) {
          df.format <- allele_format
        } else {
          df.format <- merge(df.format, allele_format, by = "Sample")
        }
      }
    }
  }
  tumor_df <- unite(df.format, "HLA Class II", 2:ncol(df.format), remove=T, sep=',')
  print(tumor_df)
  final_df <- rbind(final_df, tumor_df)
}
write_tsv(final_df, "../batch3_classII_hlatypes_unpaired.tsv")

# then in bash:
# IFS=$'\n'; for line in $(<batch3_classII_hlatypes_unpaired.tsv); do echo -e "$(echo $line | cut -f1)\t$(echo $line | cut -f2 | tr ',' '\t')"; done > batch3_classII_hlatypes_unpaired.format.tsv

# then after manual pairing:
# IFS=$'\n'; for line in $(<batch3_classII_hlatypes_paired.tsv); do echo -e "$(echo $line | cut -f1)\t$(echo $line | cut -f2- | tr '\t' ',')"; done  > batch3_classII_hlatypes_paired.format.tsv










