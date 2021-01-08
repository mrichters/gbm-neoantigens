# 8.28.19

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

get_gbm_order <- function(vector) {
  re_samples <- unique(vector[grepl("Re", vector)])
  not_re_samples <- unique(vector[!grepl("Re", vector)])
  sample_order <- c(not_re_samples, re_samples)
  tumor_order <- unique(sapply(sample_order, function(x) return(substr(x, 1, nchar(x)-2))))
  order_list <- list(sample_order, tumor_order)
  return(sample_order)
}

# set wd
setwd("~/Desktop/GBM\ Figures")

# creating manual color palette for plotting
my_pal <- c(pal_npg("nrc", alpha = 1)(10), pal_locuszoom("default")(7), pal_lancet("lanonc")(9), pal_simpsons("springfield", alpha = 1)(10))
# for y axis range!!
#+ scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))

# list all sheets in excel file
#twck_samples <- excel_sheets(path="~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx")
#twck_tumors <- unique(sapply(twck_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
#htc_samples <- excel_sheets(path="~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx")
#htc_tumors <- unique(sapply(htc_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
#variants <- "~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx"
#variants <- "~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx"

variant_path <- "~/Desktop/GBM/gbm_google_drive/merged_annotated_variants_with_primary.xlsx"
samples <- excel_sheets(path = variant_path)
tumors <- unique(sapply(samples, function(x) return(unlist(strsplit(x, '-'))[1])))

total_df <- c()
for ( tumor in tumors ) {
  tumor_df <- c()
  sample_names <- samples[grep(tumor, samples)]
  for ( name in samples[grep(tumor, samples)] ) {
    # read in correct batch variant excel file
    variants_df <- as.data.frame(read_excel(variant_path, sheet = name, col_names=T))
    # if ( grepl("Re", name) == FALSE ) {
    #   new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
    #   tumor <- unlist(str_split(name, '-'))[2]
    # } else { 
    #   new_name <- paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
    #   tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
    # } 
    # sample_names <- c(sample_names, new_name)
    # for creating master figure w/ mutational signatures: need PROPORTION of each variant type that i'm getting from the variants_df
    # for df, every consequence is a column and value is relative proportion (count / total)
    #variants_by_conseq = 0
    vbc_vector <- c()
    for ( consequence in variants_df$Consequence ) {
      if ( grepl("missense", consequence) || grepl("inframe", consequence) || grepl("frameshift", consequence) || grepl("protein_altering", consequence) && grepl("stop_gained", consequence) == FALSE) {
        if ( grepl("stop_gained", consequence) == FALSE ) {
        if ( grepl("&", consequence) ) {
          consequence <- unlist(strsplit(consequence, "&"))[1]
        }
        vbc_vector <- c(vbc_vector, consequence)
      } 
    }
  }
    sample_df <- cbind(as.data.frame(table(vbc_vector)))
    sample_df$Samples <- new_name
    sample_df$Tumor <- tumor
    tumor_df <- rbind(tumor_df, sample_df)
  }
total_df <- rbind(total_df, tumor_df)
}
colnames(total_df) <- c("Variant Consequence", "Variant Count", "Samples", "Tumor")
# for proportions
#function
#df[, -1] <- lapply( df[ , -1], function(x) x/sum(x, na.rm=TRUE) )
total_df.spread <- total_df[1:3] %>%
  spread(key = Samples, value = `Variant Count`)

total_df.spread[is.na(total_df.spread)] <- 0

total_df.spread[, -1] <- lapply(total_df.spread[ , -1], function(x) x/sum(x, na.rm=TRUE) )

total_df.tidy <- total_df.spread %>%
  gather(key = "Samples", value = "Variant Proportion", -`Variant Consequence`)

# raw counts
ggplot(total_df, aes(x = Samples, y = `Variant Count`, fill = `Variant Consequence`)) + geom_bar(stat="identity") + scale_fill_manual(values = my_pal[c(3,4,1,2,8,9)]) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14))

# proportion
ggplot(total_df.tidy, aes(x = Samples, y = `Variant Proportion`, fill = `Variant Consequence`)) + geom_bar(stat="identity") + scale_fill_manual(values = my_pal[c(3,4,1,2,8,9)]) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5), text = element_text(size=14))

### MISSENSE MUTATION SPECTRUM ###
missense_df <- c()
total_df <- c()
for ( tumor in tumors ) {
  tumor_df <- c()
  sample_names <- samples[grep(tumor, samples)]
  for ( name in sample_names ) {
    # read in correct batch variant excel file
    variants_df <- as.data.frame(read_excel(variant_path, sheet = name, col_names=T))
    missense_count = 0
    vbc_vector <- data.frame()
    for ( line in 1:nrow(variants_df) ) {
      if ( grepl("missense", variants_df$Consequence[line]) ) {
        missense_count <- missense_count + 1
        if ( grepl("&", variants_df$Consequence[line]) ) {
          consequence <- unlist(strsplit(variants_df$Consequence[line], "&"))[1]
        } else {
          consequence <- variants_df$Consequence[line]
        }
        genotype <- as.character(paste(variants_df$REF[line], variants_df$ALT[line], sep = ">"))
        if ( genotype %in% c("A>G", "C>T", "A>C", "C>G", "C>A", "A>T") ) {
          new_entry <- data.frame(consequence, genotype, name, tumor)
          vbc_vector <- rbind(vbc_vector, new_entry)
        } 
      }
    }
    if ( nrow(vbc_vector) > 0 && missense_count >= 20) {
      sample_df <- as.data.frame(table(vbc_vector))
      tumor_df <- rbind(tumor_df, sample_df)
    }
    print(missense_count)
    missense_df <- rbind(missense_df, data.frame("sample" = name, "count" = missense_count))
    missense_count = 0
  }
  total_df <- rbind(total_df, tumor_df)
}

colnames(total_df) <- c("Variant Consequence", "Variant Genotype", "Samples", "Tumor", "Frequency")

# for proportions
total_df.spread <- total_df[, c(2,3,5)] %>%
  spread(key = `Variant Genotype`, value = Frequency)

#total_df.spread$`Variant Genotype` <- sort(total_df.spread$`Variant Genotype`)
total_df.spread[is.na(total_df.spread)] <- 0
# create levels
mut_levels <- c("A>G", "C>T", "A>C", "C>G", "C>A", "A>T")
# now I need to change the level of the actual "Frequency" 
total_df.reorder <- total_df.spread[, c("Samples", mut_levels)]

total_df.reorder.switch <- gather(total_df.spread, key = "Variant Genotype", value = "Frequency", -Samples)
total_df.reorder.switch2 <- spread(total_df.reorder.switch, key = Samples, value = Frequency)
total_df.reorder.switch2[ ,-1] <- lapply(total_df.reorder.switch2[ ,-1], function(x) x/sum(x, na.rm=TRUE)*100 )
total_df.tidy <- total_df.reorder.switch2 %>%
  gather(key = "Samples", value = "Frequency", -`Variant Genotype`)
total_df.tidy.spread <- spread(total_df.tidy, key = `Variant Genotype`, value = Frequency)

total_df.final <- total_df.tidy.spread[c("Samples", mut_levels)]
total_df.final2 <- gather(total_df.final, key = `Variant Genotype`, value = Frequency, -Samples)
# change the order of the legend based on proportion
total_df.final2$`Variant Genotype` <- factor(total_df.final2$`Variant Genotype`, levels = mut_levels)
# get sample order
total_df.final2$Samples <- sapply(total_df.final2$Samples, function(x) gsub("-", "  ", x))
total_df.final2$Samples <- factor(total_df.final2$Samples, levels = get_gbm_order(unique(total_df.final2$Samples)))

mut_spec <- ggplot(total_df.final2, aes(x = Samples, y = Frequency, fill = sort(`Variant Genotype`))) + geom_bar(stat = "identity") + geom_vline(xintercept = c(4.5,7.5,10.5,13.5,16.5,19.5,22.5,25.5,28.5,31.5,33.5,35.5,39.5,41.5,42.5,45.5,48.5,51.5,53.5,56.5,59.5,62.5,65.5,67.5,69.5,73.5,76.5)) + ylab("% Total Mutations") + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8), legend.title = element_blank(),legend.position = "bottom") + scale_y_continuous(limits = c(0, 101), breaks = seq(0, 100, by = 20)) + coord_cartesian(expand = F) + scale_x_discrete(labels = c("BrMET008  1" = "1", "BrMET008  3" = "3", "BrMET008  4" = "4", "BrMET009  1" = "1", "BrMET009  3" = "3", "BrMET010  1" = "1", "BrMET010  3" = "3", "BrMET019  1" = "1", "BrMET019  3" = "3", "BrMET025  1" = "1", "BrMET025  3" = "3", "BrMET018  1" = "1", "BrMET018  3" = "3", "BrMET023  1" = "1", "BrMET023  3" = "3", "BrMET024  1" = "1", "BrMET024  3" = "3", "BrMET027  1" = "1", "BrMET027  3" = "3", "BrMET028  1" = "1", "BrMET028  3" = "3", "BrMET058  4" = "4", "GBM030  2" = "2","GBM032  1" = "1", "GBM032  3" = "3", "GBM032  4" = "4", "GBM051  1" = "1", "GBM051  3" = "3", "GBM052  3" = "3", "GBM055  1" = "1", "GBM055  3" = "3", "GBM056  1" = "1", "GBM056  3" = "3", "GBM059  1" = "1", "GBM059  3" = "3", "GBM062  1" = "1", "GBM062  3" = "3", "GBM063  1" = "1", "GBM063  3" = "3", "GBM064  1" = "1", "GBM064  3" = "3", "GBM069  1" = "1", "GBM069  3" = "3", "GBM070  1" = "1", "GBM070  3" = "3", "GBM074  1" = "1", "GBM074  3" = "3", "GBM079  1" = "1", "GBM079  3" = "3", "GBM083  1" = "1", "GBM083  3" = "3", "GBM018.Re  3" = "3", "GBM031.Re  1" = "1", "GBM031.Re  3" = "3", "GBM031.Re  4" = "4", "GBM047.Re  1" = "1", "GBM047.Re  3" = "3", "GBM065.Re  1" = "1", "GBM065.Re  3" = "3", "GBM065.Re  4" = "4"))



ggsave(mut_spec, file="mutation_spectrum.pdf", width = 7, height = 4)
 