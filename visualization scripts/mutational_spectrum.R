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
# set wd
setwd("~/Desktop/GBM_plots")

# creating manual color palette for plotting
my_pal <- c(pal_npg("nrc", alpha = 1)(10), pal_locuszoom("default")(7), pal_lancet("lanonc")(9), pal_simpsons("springfield", alpha = 1)(10))
# for y axis range!!
#+ scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))

# list 3 gbm batches for analysis
gbm_batches <- c("Batch1", "Batch2", "Batch3")

# list all sheets in excel file
#twck_samples <- excel_sheets(path="~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx")
#twck_tumors <- unique(sapply(twck_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
#htc_samples <- excel_sheets(path="~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx")
#htc_tumors <- unique(sapply(htc_samples, function(v) return(unlist(strsplit(v, '-'))[2])))
#variants <- "~/Desktop/GBM/gbm_google_drive/HTC_annotated_variants.xlsx"
#variants <- "~/Desktop/GBM/gbm_google_drive/TWCK_annotated_variants.xlsx"

total_df <- c()
for ( tumor in twck_tumors ) {
  tumor_df <- c()
  sample_names <- c()
  for ( name in twck_samples[grep(tumor, twck_samples)] ) {
    variants_df <- as.data.frame(read_excel(variants, sheet = name))
    if ( grepl("Re", name) == FALSE ) {
      new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
      tumor <- unlist(str_split(name, '-'))[2]
    } else { 
      new_name <- paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
      tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
    } 
    sample_names <- c(sample_names, new_name)
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
for ( batch in gbm_batches ) {
  variant_path <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_annotated_variants.xlsx", sep="")
  samples <- excel_sheets(path = variant_path)
  tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[2])))
  for ( tumor in tumors ) {
    tumor_df <- c()
    for ( name in samples[grep(tumor, samples)] ) {
      # read in correct batch variant excel file
      variants_df <- as.data.frame(read_excel(variant_path, sheet = name))
      # create standardized shorter sample names
      if ( grepl("Re", name) == TRUE ) {
        if ( grepl("Primary", name) == TRUE ) {
          new_name <- "Re.GBM065-P"
          tumor <- "GBM065"
        } else {
          new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
          tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
        }
      } else if ( batch == "Batch3") { 
        new_name <- mgsub(name, c("19-", "Tumor", "_"), c("", "", "-"))
        tumor <- unlist(strsplit(new_name, '-'))[1]
      } else if ( grepl("Re", name) == FALSE ) {
        new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
        tumor <- unlist(str_split(name, '-'))[2]
      } 
      #sample_names <- c(sample_names, new_name)
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
            new_entry <- data.frame(consequence, genotype, new_name, tumor)
            vbc_vector <- rbind(vbc_vector, new_entry)
          } 
        }
      }
      if ( nrow(vbc_vector) > 0 && missense_count >= 20) {
        sample_df <- as.data.frame(table(vbc_vector))
        tumor_df <- rbind(tumor_df, sample_df)
      }
      print(missense_count)
      missense_df <- rbind(missense_df, data.frame("sample" = new_name, "count" = missense_count))
      missense_count = 0
    }
    total_df <- rbind(total_df, tumor_df)
  }
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
total_df.final2$`Variant Genotype` <- as.factor(total_df.final2$`Variant Genotype`)
levels(total_df.final2$`Variant Genotype`) <- mut_levels

mut_spec <- ggplot(total_df.final2, aes(x = Samples, y = Frequency, fill = sort(`Variant Genotype`))) + geom_bar(stat = "identity") + ylab("% Total Mutations") + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=12), legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) + ggtitle("Mutation Spectrum") + scale_y_continuous(limits = c(0, 101), breaks = seq(0, 100, by = 20))



ggsave(mut_spec, file="mutation_spectrum.pdf", width = 12, height = 5)
 