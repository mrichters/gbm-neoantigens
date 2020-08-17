library(tidyverse)
library(readxl)
library(openxlsx)

setwd("~/Desktop/GBM_Project/R scripts")

----------------------------
## 5.15.19 ##
# making new variant plots w/ twck and updated htc variant counts
new_var_df <- as.data.frame(read_excel('../variant_counts.xlsx', sheet = 'reformatted variants'))
new_var_df <- as.data.frame(read_excel('../formatted_variant_counts.xlsx'))

tumor_vector <- c()
for (sample in new_var_df$Sample) {
  tumor <- unlist(strsplit(sample, "-"))[1]
  tumor_vector <- c(tumor_vector, tumor)
}
new_var_df$Tumor <- tumor_vector
bm_re_samples <- new_var_df[c(1:29), ]
gbm_samples <- new_var_df[c(30:59), ]
#barplot
ggplot(gbm_samples, aes(x=Sample, y=`Number of Variants`, fill=Tumor)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("Total Variants") #+ scale_fill_brewer(palette="Dark2")
###
#all formatting commands for variant_counts.xlsx --> I manually reformatted this file into formatted_variant_counts.xlsx for plotting
tumor_vector <- c()
sample_vector <- c()
for ( sample in new_var_df$Sample) {
  sample <- gsub("H_TC", "HTC", sample)
  new_sample <- gsub("_", "-", sample)
  split_sample <- unlist(strsplit(new_sample, "-"))
  tumor <- split_sample[2]
  sample_num <- split_sample[5]
  tumor_vector <- c(tumor_vector, tumor)
  sample_name <- paste(tumor, sample_num, sep='-')
  sample_vector <- c(sample_vector, sample_name)
} 
new_var_df$Tumor <- tumor_vector
new_var_df$`Sample Name` <- sample_vector
twck <- new_var_df[23:61,]
recurrence <- twck[grep(".Re", twck$`Sample Name`),]
no_recurrence <- twck[-grep(".Re", twck$`Sample Name`),]
brain_mets <- twck[grep("BrMET", twck$`Sample Name`),]
no_brain_mets <- twck[-grep("BrMET", twck$`Sample Name`),]
gbm_only <- no_recurrence[-grep("BrMET", twck$`Sample Name`),]
ordered_new <- rbind(brain_mets, recurrence, gbm_only)
###

#barplot
ggplot(new_var_df, aes(x=Sample, y=`Number of Variants`, fill=Tumor)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("Total Variants") + ylim(0,300)#+ scale_fill_brewer(palette="Dark2")

------------------------------

## NUMBER OF EXPRESSED/TOTAL VARIANTS AND NEOANTIGENS PER SAMPLE ##

# for variant - use 'gbm_google_drive/GBM_variants_per_sample.xlsx'
# for neoantigens class I - use 'gbm_google_drive/GBM_pvactools_clinicalclassI_11.20.18.xlsx'
# for neoantigens class II - use 'gbm_google_drive/GBM_pvactools_clinicalclassII_condensed_11.28.18.xlsx'

# create sample names vectors!
HTC_sample_names <- excel_sheets(path='../gbm_google_drive/GBM_HTC_annotated_variants.xlsx')
HTC_tumors <- sapply(HTC_sample_names, function(v) return(unlist(strsplit(v, '-'))[2]))
TWCK_sample_names <- excel_sheets(path='../gbm_google_drive/GBM_TWCK_annotated_variants.xlsx')
TWCK_tumors <- sapply(TWCK_sample_names, function(v) return(unlist(strsplit(v, '-'))[2]))


## EXPRESSED ##

# for loop for variants (gene expression info is harder to parse)
df.variants <- data.frame()
df.totalgx <- data.frame()
for (i in HTC_sample_names) {
  df <- as.data.frame(read_excel('../gbm_google_drive/GBM_HTC_annotated_variants.xlsx', sheet=i))
  df.subset <- df[ ,c(1,2,15)] #for variants
  colnames(df.subset) <- c("CHROM", "POS", "Gene_Expression")
  df.subset[is.na(df.subset)] <- 0
  df.gx <- data.frame()
  df.sample <- data.frame()
  for ( row in 1:nrow(df.subset) ) {
    gx <- as.character(df.subset[row, "Gene_Expression"])
    chrom <- df.subset[row, "CHROM"]
    pos <- df.subset[row, "POS"]
    if ( gx == "0" ) {
      gx.exp <- as.numeric(gx)
      gx.gene <- NA
    } else if ( grepl(',', gx) == TRUE ) {
      gx <- strsplit(gx, '\\,')
      gx <- unlist(gx)[1]
      gx <- strsplit(gx, '\\|')
      gx.exp <- as.numeric(unlist(gx)[2])
      gx.gene <- unlist(gx)[1]
    } else if ( grepl(',', gx) == FALSE ) {
      gx <- strsplit(gx, '\\|')
      gx.exp <- as.numeric(unlist(gx)[2])
      gx.gene <- unlist(gx)[1]
      print(gx.exp)
    }
    if ( gx.exp != 0) {
      df_append <- data.frame(i, chrom, pos, gx.gene, gx.exp)
      df.gx <- rbind(df.gx, df_append) 
    }
  }
  Variants <- nrow(df.gx)
  Samples <- i
  df.sample <- data.frame(Samples, Variants)
  df.variants <- rbind(df.variants, df.sample)
  df.totalgx <- rbind(df.totalgx, df.gx)
} 

# for neoantigens
df.na <- data.frame()
df.totalna_I <- data.frame()
for (i in HTC_sample_names) {
  df <- as.data.frame(read_excel('../gbm_google_drive/GBM_pvactools_clinicalclassI_condensed_11.20.18.xlsx', sheet=i))
  df.subset <- df[ ,c(7,15)] #for neoantigens
  colnames(df.subset) <- c("Neoantigen", "Gene_Expression")
  df.subset[is.na(df.subset)] <- 0
  df.gx <- data.frame()
  df.sample <- data.frame()
  if ( nrow(df.subset) > 0 ) {
    for ( row in 1:nrow(df.subset) ) {
      gx <- df.subset[row, "Gene_Expression"]
      neoantigen <- df.subset[row, "Neoantigen"]
      if ( gx > 0 ) {
        df_append <- data.frame(i, neoantigen, gx)
        df.gx <- rbind(df.gx, df_append) 
      }
    }
  } else {
    df.gx <- data.frame()
}
  Neoantigens <- nrow(df.gx)
  Samples <- i
  df.sample <- data.frame(Samples, Neoantigens)
  df.totalna_I <- rbind(df.totalna_I, df.sample)
} 
df.totalna_I$Tumor <- c(rep("BrMET008", times=4), rep("BrMET009", times=3), rep("GBM030", times=2), rep("GBM032", times=3), rep("GBM047", times=2), rep("GBM051", times=3), rep("GBM052", times=2))
df.totalna_I$Samples <- expressed_sample_names
exp_df.combined.noexp <- df.combined.noexp[c(1,2,3,4,5,6,7,8,9,10,11,12,15,16,17,18,19,21,22), ]
colnames(exp_df.combined.noexp)[1] <- "Samples"
exp_df.combined <- merge(df.variants, df.totalna_I, df.totalna_II)

na_I_combined <- merge(exp_df.combined.noexp, df.totalna_I, by="Samples")

----------------------------------------
# 5.20.19 #
## NOT EXPRESSED ##
# HTC #
df.htc <- data.frame()
for ( i in HTC_sample_names ) {
  df_var <- as.data.frame(read_excel('../gbm_google_drive/GBM_HTC_annotated_variants.xlsx', sheet=i))
  df_fil <- as.data.frame(read_excel('../gbm_google_drive/GBM_variants_per_sample_v2.xlsx', sheet=i))
  #df_na_I <- as.data.frame(read_excel('../gbm_google_drive/GBM_pvactools_clinicalclassI_11.20.18.xlsx', sheet=i))
  #df_na_II <- as.data.frame(read_excel('../gbm_google_drive/GBM_pvactools_clinicalclassII_condensed_11.28.18.xlsx', sheet=i))
  Variants_2019 <- nrow(df_var)
  Variants_2018 <- nrow(df_fil)
  #ClassI_Neoantigens <- nrow(df_na_I)
  #ClassII_Neoantigens <- nrow(df_na_II)
  Sample <- i
  Tumor <- unlist(strsplit(i, '-'))[2]
  Short_sample <- paste(Tumor, substr(Sample, nchar(Sample), nchar(Sample)), sep="-")
  #df_append <- data.frame(i, ClassI_Neoantigens, ClassII_Neoantigens, Filtered_Variants, Variants)
  #df.combined.noexp <- rbind(df.combined.noexp, df_append)
  df_append <- data.frame(Sample, Short_sample, Tumor, Variants_2018, Variants_2019)
  df.htc <- rbind(df.htc, df_append)
}
df.htc.clean <- gather(df.htc, Analysis, Variants, 4:5)

# side by side plots
ggplot(df.htc.clean, aes(x = Short_sample, y = Variants, fill = Tumor)) + geom_bar(stat='identity') + xlab("Samples") + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=18), strip.text.x = element_text(size = 17)) + facet_grid(. ~ Analysis)

# TWCK #
TWCK_sample_names <- excel_sheets(path='../gbm_google_drive/twck_pvacseq_2019/TWCK_classI_pvacseq_results.xlsx')
df.twck <- data.frame()
for ( i in TWCK_sample_names ) {
  #df_var <- as.data.frame(read_excel('../gbm_google_drive/GBM_TWCK_annotated_variants.xlsx', sheet=i))
  df_var <- as.data.frame(read_excel('../gbm_google_drive/twck_pvacseq_2019/TWCK_classI_pvacseq_results.xlsx', sheet=i))
  Neoantigens <- nrow(df_var)
  Sample <- i
  Tumor <- unlist(strsplit(i, '-'))[2]
  if ( grepl("Re", Sample) ) {
    Short_sample <- paste(Tumor, "Re", substr(Sample, nchar(Sample), nchar(Sample)), sep="-")
  } else {
    Short_sample <- paste(Tumor, substr(Sample, nchar(Sample), nchar(Sample)), sep="-")
  } 
  df_append <- data.frame(Sample, Short_sample, Tumor, Neoantigens)
  df.twck <- rbind(df.twck, df_append)
}

# HTC #
HTC_sample_names <- excel_sheets(path='../gbm_google_drive/htc_pvacseq_2019/HTC_classII_pvacseq_results.xlsx')
df.htc.II <- data.frame()
for ( i in HTC_sample_names ) {
  #df_var <- as.data.frame(read_excel('../gbm_google_drive/GBM_TWCK_annotated_variants.xlsx', sheet=i))
  df_var <- as.data.frame(read_excel('../gbm_google_drive/htc_pvacseq_2019/HTC_classII_pvacseq_results.xlsx', sheet=i))
  Neoantigens <- nrow(df_var)
  Sample <- i
  Tumor <- unlist(strsplit(i, '-'))[2]
  if ( grepl("Re", Sample) ) {
    Short_sample <- paste(Tumor, "Re", substr(Sample, nchar(Sample), nchar(Sample)), sep="-")
  } else {
    Short_sample <- paste(Tumor, substr(Sample, nchar(Sample), nchar(Sample)), sep="-")
  } 
  df_append <- data.frame(Sample, Short_sample, Tumor, Neoantigens)
  df.htc.II <- rbind(df.htc.II, df_append)
}
# write all dfs to csv and format for final plots
write_csv(df.htc, "df_htc.csv")

#formatted info - formatted_variant_counts.xlsx
neo_df_new <- as.data.frame(read_excel('../formatted_variant_counts.xlsx'))
neo_df.tidy <- neo_df %>% gather(Variants, `Class I`, `Class II`, key = "Type", value = "Number")

bm_re <- neo_df[c(1:29),]
gbm_primary <- neo_df[c(30:59),]

neo_only_df <- gbm_primary[, c(1,2,4,5)]
neo_only_df.tidy <- neo_only_df %>% gather(`Class I`, `Class II`, key = "Neoantigens", value = "Number")

neo_df_new$`Class I` <- as.numeric(neo_df_new$`Class I`)
neo_df_new$`Class II` <- as.numeric(neo_df_new$`Class II`)
neo_df_new$`Total Neoantigens` <- neo_df_new$`Class I` + neo_df_new$`Class II`

#dot plot
ggplot(neo_df_new, aes(x=`Total Neoantigens`, y=Variants, color=Tumor)) + geom_point(size=4) + theme(text = element_text(size=16)) + xlab("Total Neoantigens") + ylim(0,400) + xlim(0,300)

# bar plot
ggplot(neo_only_df.tidy, aes(x=Sample, y=Number, fill=Neoantigens)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + xlab("Samples") + ylab("Count")


  
------------------------------------------

df.combined.noexp$Tumor <- c(rep("BrMET008", times=4), rep("BrMET009", times=3), rep("GBM030", times=2), rep("GBM032", times=3), rep("GBM047", times=2), rep("GBM051", times=3), rep("GBM052", times=2))
df.combined.noexp$Samples <- expressed_sample_names


colnames(df.variants) <- c("Samples", "Count")
df.variants$Tumor <- expressed_tumor_names
df.variants$Type <- "Variants"
df.combined.exp <- rbind(df.classI_na, df.classII_na, df.variants)
tumor_names <- c(rep("BrMET008", 4), rep("BrMET009", 3), rep("GBM030", 2), rep("GBM032", 4), rep("GBM047", 3), rep("GBM051", 3), rep("GBM052", 3))
expressed_tumor_names <- c(rep("BrMET008", 4), rep("BrMET009", 3), rep("GBM030", 2), rep("GBM032", 3), rep("GBM047", 2), rep("GBM051", 3), rep("GBM052", 2))
df.total$Tumor <- expressed_tumor_names
#  Variants <- as.numeric(nrow(df))
#  Samples <- i
#  df <- data.frame(Samples, Variants)
#  final_df_var_fil <- rbind(final_df_var_fil, df) 
#}

sample_names <- c(paste("BrMET008-", 1:4, sep = ""), paste("BrMET009-", 1:3, sep = ""), paste("GBM030-", 1:2, sep = ""), paste("GBM032-", 1:4, sep = ""), paste("GBM047-", 1:3, sep = ""), paste("GBM051-", 1:3, sep = ""), paste("GBM052-", 1:3, sep = ""))

expressed_sample_names <- c(paste("BrMET008-", 1:4, sep = ""), paste("BrMET009-", 1:3, sep = ""), paste("GBM030-", 1:2, sep = ""), paste("GBM032-", 1:3, sep = ""), paste("GBM047-", 2:3, sep = ""), paste("GBM051-", 1:3, sep = ""), paste("GBM052-", 2:3, sep = ""))

df.combined.noexp$Samples <- sample_names
df.combined.noexp$Tumor <- tumor_names
# bar plot showing number of variants per sample
display.brewer.all()

df.combined$Samples <- rep(expressed_sample_names, 3)
df.combined$Full_Samples <- rep(sample_names, 3)
df.combined$Exp_Type <- c(rep("Expressed Class I Neoantigens", 19), rep("Expressed Class II Neoantigens", 19), rep("Expressed Variants", 19))

# barplots
ggplot(new_var_df[1:22, c(1,3)], aes(x=Sample, y=ClassI_Neoantigens, fill=Tumor)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("Class I Neoantigens") + scale_fill_brewer(palette="Dark2") + ylim(0,100)

# boxplot
ggplot(df.combined, aes(x=Tumor, y=Count, color=Tumor)) + geom_boxplot() + theme(text = element_text(size=18)) + guides(color=FALSE) + ylab("Class II Neoantigens") + scale_fill_brewer(palette="Dark2")

# side by side plots
ggplot(df.combined, aes(x = Tumor, y = Count, col = Tumor)) + geom_boxplot(size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=18), strip.text.x = element_text(size = 17)) + facet_grid(. ~ Exp_Type) + guides(col=FALSE) + scale_fill_brewer(palette="Dark2")
#ggplot(df.combined, aes(x = Samples, y = Count, col = Tumor)) + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=18), strip.text.x = element_text(size = 17)) + facet_grid(. ~ Exp_Type) + guides(col=FALSE) + scale_fill_brewer(palette="Dark2")




# combine the data.frames
combined_df <- merge(final_df_var, final_df_var_fil, by="Samples")
#combined_df$Total_Neoantigens <- rowSums(combined_df[,c(2,3)])
combined_df$ClassI <- final_df_classI$Variants
combined_df$ClassII <- final_df_classII$Variants
colnames(combined_df) <- c("Samples", "Total_Variants", "Filtered_Variants","ClassI_Neoantigens", "ClassII_Neoantigens")
# add color dimension
combined_df$Tumor <- c(rep("BrMET008", times=4), rep("BrMET009", times=3), rep("GBM030", times=2), rep("GBM032", times=4), rep("GBM047", times=3), rep("GBM051", times=3), rep("GBM052", times=3))
#combined_df$Tumor_xaxis <- c(rep("008", times=4), rep("009", times=3), rep("030", times=2), rep("032", times=4), rep("047", times=3), rep("051", times=3), rep("052", times=3))
# tidy the data for side by side plots
combined_df.tidy <- gather(combined_df, Type, Count, Total_Variants, Filtered_Variants, ClassI_Neoantigens, ClassII_Neoantigens, -Samples, -Tumor)
combined_df.tidy$Type <- as.factor(combined_df.tidy$Type)

# plot the data! - scatter
ggplot(combined_df, aes(x=Total_Neoantigens, y=Variants, color = Tumor)) + geom_point(size = 4.5, alpha = 0.8) + theme(text = element_text(size=15)) 
# boxplot
ggplot(df.totalna_II, aes(x=Tumor, y=ClassII, color=Tumor)) + geom_boxplot() + theme(text = element_text(size=18)) + guides(color=FALSE) + ylab("Class II Neoantigens") + ylim(0,200)
# side by side plots
ggplot(combined_df.tidy, aes(x = Tumor, y = Count, col = Tumor)) + geom_boxplot(size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=18), strip.text.x = element_text(size = 17)) + facet_grid(. ~ Type) + guides(col=FALSE)





---------------------------------------
### create scatter plot for how many variants per sample create neoantigens ###
setwd("~/Desktop/GBM")
gbm_batches <- c("Batch1", "Batch2", "Batch3")

percentage_df <- data.frame()
for ( batch in gbm_batches ) {
    variant_file <- paste("gbm_google_drive/", batch, "_annotated_variants.xlsx", sep="")
    naI_file <- paste("/Users/mrichters/Desktop/GBM/gbm_google_drive/", batch, "_pvacseq/clinical/", batch, "_clinical_classI_pvacseq_ranked.xlsx", sep="")
    naII_file <- paste("/Users/mrichters/Desktop/GBM/gbm_google_drive/", batch, "_pvacseq/clinical/", batch, "_clinical_classII_pvacseq_ranked.xlsx", sep="")
    sample_names <- excel_sheets(variant_file)
    for (i in sample_names) {
        df_var <- as.data.frame(read_excel(variant_file, sheet=i))
        df_na_I <- as.data.frame(read_excel(naI_file, sheet=i))
        df_na_II <- as.data.frame(read_excel(naII_file, sheet=i))
        variant_count <- nrow(df_var)
        if ( variant_count == 0 || dim(df_na_I)[1] == 0 && dim(df_na_II)[1] == 0 ) {
            percentage_line <- data.frame("Sample" = i, "Total_Variants" = nrow(df_var), "Total_Neoantigens" = nrow(df_na_I)+nrow(df_na_II), "Total_Neoantigen-forming_Variants" = 0, "Proportion" = 0)
            percentage_df <- rbind(percentage_df, percentage_line)
        } else {
            variant_count <- nrow(df_var)
            neoantigens <- rbind(df_na_I[, c("HGVSc", "HLA Allele", "MT Epitope Seq")], df_na_II[, c("HGVSc", "HLA Allele", "MT Epitope Seq")])
            neoantigens$Sample <- i
            #na_df <- rbind(na_df, neoantigens) 
            neoantigen_count <- nrow(neoantigens)
            na_variants <- unique(sort(as.character(neoantigens[, 1])))
            na_variants_count <- length(na_variants)
            percentage_line <- data.frame("Sample" = i, "Total_Variants" = variant_count, "Total_Neoantigens" = neoantigen_count, "Total_Neoantigen-forming_Variants" = na_variants_count, "Proportion" = na_variants_count/variant_count)
            percentage_df <- rbind(percentage_df, percentage_line)
        } 
    }
}
# write df to excel file
write.xlsx(percentage_df, "neoantigen_forming_variants_proportion.xlsx")

# plot
ggplot(percentage_df, aes(x=Total_Variants, y=Total_Neoantigens)) + geom_point(size = 2, color="grey20") + theme(text = element_text(size=14), panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + geom_smooth(method=lm) + ylab("Total Neoantigens") + xlab("Total Variants")

## FIND R SQUARED VALUE OF VARIANTS / NEOANTIGENS CORRELATION
percentage_df.lm = lm(Total_Variants ~ Total_Neoantigen.forming_Variants, data=percentage_df)
summary(percentage_df.lm)$adj.r.squared
summary(percentage_df.lm) #<-- tells us p-value and whether correlation is significant

purity_check <- as.data.frame(read_excel("neoantigen_forming_variants_proportion.xlsx"))
purity_check$`Tumor Purity` <- as.numeric(purity_check$`Tumor Purity`)
ggplot(purity_check, aes(x=Proportion, y=`Tumor Purity`)) + geom_point(size = 3) + theme(text = element_text(size=16)) + geom_smooth(method=lm)

brmets <- filter(purity_check, str_detect(Sample, "BrMET"))
gbms <- filter(purity_check, str_detect(Sample, "GBM"))

x <- brmets
print(c(min(x$Total_Variants), max(x$Total_Variants), median(x$Total_Variants)))

### 1.15.2020 ###
my_pal <- c(pal_d3("category10", alpha = 1)(10), pal_aaas("default", alpha = 1)(7), pal_d3("category20c", alpha = 1)(10))

purity_check$Tumor <- sapply(purity_check$Sample, function(v) return(unlist(strsplit(v, '-'))[2])) 

tumor_type_vector <- c()
for ( row in 1:nrow(purity_check) ) {
  tumor <- purity_check[row, "Tumor"]
  if ( grepl("BrMET", tumor) ) {
    tumor_type <- "BrMET"
  } else {
    tumor_type <- "GBM"
  }
  tumor_type_vector <- c(tumor_type_vector, tumor_type)
}
purity_check$Tumor_Type <- tumor_type_vector

ggplot(purity_check, aes(x=Tumor_Type, y=Proportion)) + geom_boxplot(notch=FALSE, outlier.shape = NA, width=0.3) + geom_jitter(position=position_jitter(0.15),aes(color=Tumor_Type), size=2) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=14), legend.position = "none", panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + ylab("Variants resulting in neoantigens") + xlab("Tumor Type") + scale_color_manual(values = c("indianred4", "grey30"))

  