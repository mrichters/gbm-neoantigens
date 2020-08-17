# **** new seq2hla code at bottom
# import packages
library(readxl)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggsci)

excel_file <- '../hla_abundance_files/GBM_hla_expression.xlsx'
sheet_names <- excel_sheets(excel_file)
tumor_names <- c("BrMET008", "BrMET009", "GBM030", "GBM032", "GBM047", "GBM051", "GBM052")

### CLASS I HLA ALLELES ###

for ( tumor in tumor_names ) {
  combined_df <- data.frame()
  for ( sheet in excel_sheets(excel_file) ) {
    if ( grepl(tumor,sheet) == TRUE ) {
      df <- as.data.frame(read_excel(excel_file, sheet=sheet))
      tumor_num <- paste("-Tumor-", substr(sheet, nchar(sheet), nchar(sheet)), sep="") 
      df$Sample <- paste(tumor, tumor_num, sep="")
      if ( tumor == "GBM030" || tumor == "GBM032") {
        classI_df <- df[1:5, ]
      } else {
        classI_df <- df[1:6, ]
      }
      classI_df <- classI_df[order(classI_df$target_id),]
      alleles_vector <- character()
      generic_alleles_vector <- character()
      orig_generic_alleles_vector <- character()
      hla_genes_vector <- character()
      counter <- 1
      for ( item in classI_df$target_id ) {
        new_item <- sub("_", "-", substr(item, 1, 7))
        new_generic_item <- substr(item, 1, 1)
        if ( is.element(new_generic_item, orig_generic_alleles_vector) ) {
          counter <- counter + 1
          final_generic_item <- paste(new_generic_item, counter, sep="")
        } else {
          counter <- 1
          final_generic_item <- paste(new_generic_item, counter, sep="")
          orig_generic_alleles_vector <- c(orig_generic_alleles_vector, new_generic_item)
        }
        alleles_vector <- c(alleles_vector, new_item)
        generic_alleles_vector <- c(generic_alleles_vector, final_generic_item)
        hla_genes_vector <- c(hla_genes_vector, new_generic_item)
      }
      classI_df$`HLA Alleles` <- alleles_vector
      classI_df$generic_hla_alleles <- generic_alleles_vector
      classI_df$Tumor <- tumor
      classI_df$`HLA Genes` <- hla_genes_vector
      combined_df <- rbind(combined_df, classI_df)
    }
    assign(paste(tumor, '_classI', sep=""), combined_df)
  }
}

# combine all data
facet_wrap_plot <- rbind(BrMET008_classI, BrMET009_classI, GBM030_classI, GBM032_classI, GBM047_classI, GBM051_classI, GBM052_classI)

# subset the data
#facet_A <- facet_wrap_plot[facet_wrap_plot$generic_hla_alleles == 'A1',]
facet_A <- facet_wrap_plot[grep("^A", facet_wrap_plot$`HLA Alleles`), ]
facet_A$`HLA Genes` <- 'HLA-A'
  
# facet wrap barplot
ggplot(facet_wrap_plot, aes(x=generic_hla_alleles, y=tpm, fill=Sample)) + geom_bar(stat='identity', position='dodge') + facet_wrap(. ~ Tumor) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("TPM") + xlab("Class I HLA Alleles") 

# individual plots
ggplot(GBM052_classI, aes(x=HLA_Alleles, y=tpm, fill=Sample)) + geom_bar(stat='identity', position='dodge') + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("TPM") + xlab("Class I HLA Alleles") + scale_fill_brewer(palette="Dark2")

# save plot data
write.csv(facet_wrap_plot, 'classI_hla_expression_plot_data.csv')

# 3.28.19 - separate allele plots
# violin plot
ggplot(facet_C, aes(x=`HLA Genes`, y=tpm)) + geom_violin() + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("TPM") + geom_jitter(data = facet_C, aes(col = Tumor, fill = Tumor), size = 2.5) + scale_shape_manual(values=c(15,16,17,3,4,23,25,7,8,9)) + theme(legend.text = element_text(size = 9)) +  theme(legend.title = element_text(size = 10, face='bold')) #+ theme(legend.position="none")
#aes(col = Tumor, shape = `HLA Alleles`, fill = Tumor)

# facet grid plot
#ggplot(facet_C, aes(x=generic_hla_alleles, y=tpm, fill=Sample)) + geom_bar(stat='identity', position='dodge') + facet_grid(. ~ Tumor) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("TPM") + xlab("HLA-A Alleles") 


### CLASS II HLA ALLELES ###

for ( tumor in tumor_names ) {
  combined_df <- data.frame()
  for ( sheet in excel_sheets(excel_file) ) {
    if ( grepl(tumor,sheet) == TRUE ) {
      df <- as.data.frame(read_excel(excel_file, sheet=sheet))
      tumor_num <- paste("-Tumor-", substr(sheet, nchar(sheet), nchar(sheet)), sep="") 
      df$Sample <- paste(tumor, tumor_num, sep="")
      classII_df <- df[7:nrow(df), ]
      classII_df <- classII_df[order(classII_df$target_id),]
      alleles_vector <- character()
      generic_alleles_vector <- character()
      orig_generic_alleles_vector <- character()
      hla_genes_vector <- character()
      hla_gene_family_vector <- character()
      counter <- 1
      for ( item in classII_df$target_id ) {
        new_item <- sub(":", "-", substr(item, 1, 10))
        new_generic_item <- substr(item, 1, 4)
        gene_item <- substr(item, 1, 2)
        if ( is.element(new_generic_item, orig_generic_alleles_vector) ) {
          counter <- counter + 1
          final_generic_item <- paste(new_generic_item, counter, sep="-")
          #orig_generic_alleles_vector <- c(orig_generic_alleles_vector, new_generic_item)
        } else {
          counter <- 1
          final_generic_item <- paste(new_generic_item, counter, sep="-")
          orig_generic_alleles_vector <- c(orig_generic_alleles_vector, new_generic_item)
        }
        alleles_vector <- c(alleles_vector, new_item)
        generic_alleles_vector <- c(generic_alleles_vector, final_generic_item)
        hla_genes_vector <- c(hla_genes_vector, new_generic_item)
        hla_gene_family_vector <- c(hla_gene_family_vector, gene_item)
      }
      classII_df$`HLA Alleles` <- alleles_vector
      classII_df$generic_hla_alleles <- generic_alleles_vector
      classII_df$Tumor <- tumor
      classII_df$`HLA Genes` <- hla_genes_vector
      classII_df$`HLA Gene Class` <- hla_gene_family_vector
      combined_df <- rbind(combined_df, classII_df)
    }
    assign(paste(tumor, "_classII", sep=""), combined_df)
  }
}
# combine all data
facet_wrap_plot_classII <- rbind(BrMET008_classII, BrMET009_classII, GBM030_classII, GBM032_classII, GBM047_classII, GBM051_classII, GBM052_classII)

# subset the class II alleles to make facet wrap plot more readable
DP <- facet_wrap_plot_classII[grep("^DPB1", facet_wrap_plot_classII$`target_id`), ]

# facet wrap barplot
ggplot(DR, aes(x=generic_hla_alleles, y=tpm, fill=Sample)) + geom_bar(stat='identity', position='dodge') + facet_wrap(. ~ Tumor) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("TPM") + xlab("HLA-DR Class II Alleles")

# 4.1.19
# violin plot
ggplot(DP, aes(x=`HLA Genes`, y=tpm)) + geom_violin() + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("TPM") + geom_jitter(data = DP, aes(col = Tumor, fill = Tumor), size = 2.5) + scale_shape_manual(values=c(15,16,17,3,4,23,25,7,8,9)) + theme(legend.text = element_text(size = 9)) +  theme(legend.title = element_text(size = 10, face='bold')) #+ theme(legend.position="none")

# save plot data
write.csv(facet_wrap_plot, 'classII_hla_expression_plot_data.csv')

#######################################

# 8.7.19
# using seq2hla tsv combined output file

# set wd
setwd("~/Desktop/GBM_plots/most_recent_plots/seq2hla_hla_expression_jitter_plots/")

# load in file
df <- read.table("~/Desktop/GBM/seq2hla_expression/gbm_hla_expression_seq2hla.tsv", header = TRUE)
df2 <- read.table("~/Desktop/GBM/seq2hla_expression/twck_hla_expression_seq2hla.tsv", header = TRUE)
# shorten sample names / create tumor names for plotting
Samples <- c()
Tumors <- c()

for ( name in df[, 1] ) {
  if ( grepl("Re", name) == FALSE ) {
    new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
    tumor <- unlist(str_split(name, '-'))[2]
  } else {
    new_name <- paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
    tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")  
  } 
  Samples <- c(Samples, new_name)
  Tumors <- c(Tumors, tumor)
}
df$Samples <- Samples
df$Tumor <- Tumors
df$Batch <- c(rep("Batch 1", 19), rep("Batch 2", 33))
# sort the data by sample type (brmet, gbm, re.gbm)
df.order <- df[order(df$Samples),]

for ( line in 1:nrow(df.order) ) {
  if ( grepl("BrMET", df.order$Tumor[line]) == TRUE ) {
    df.order$`Tumor Type`[line] <- "Brain Metastasis"
  } else if ( grepl("Re", df.order$Tumor[line]) == TRUE ) {
    df.order$`Tumor Type`[line] <- "Recurrent GBM"
  } else {
    df.order$`Tumor Type`[line] <- "Primary GBM"
  }
}

# tidy the data
df.tidy <- gather(df.order, key = HLA_Allele, value = RPKM, -Sample, -Samples, -Tumor, -`Tumor Type`, -Batch)
df.tidy$`HLA Class` <- c(rep("Class I", 156), rep("Class II", 312))
# class I subset
df.subset.I <- df.tidy[1:156, ]
# class II subset
df.subset.II <- df.tidy[157:468, ]

# plotting data in df.order

# facet wrap barplot
ggplot(df.subset.II, aes(x=Samples, y=RPKM, fill=HLA_Allele)) + geom_bar(stat='identity', position='dodge') + facet_grid(. ~ HLA_Allele) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + theme(axis.text.x=element_blank())

# jitter + boxplot colored by HLA Allele
my_pal <- c(pal_npg("nrc", alpha = 1)(10), pal_simpsons("springfield", alpha = 1)(10))

pdf("combined_classII.pdf", width = 8, height = 8)
ggplot(df.subset.I, aes(x=HLA_Allele, y=RPKM)) + geom_boxplot(notch=FALSE, outlier.shape = NA, width=0.4) + geom_jitter(position=position_jitter(0.15),aes(color=Tumor), size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + ylab("RPKM") + scale_color_manual(values = my_pal) + xlab("Class I HLA Alleles")
plot(my_plot)
dev.off()

# same but colored by tumor type
ggplot(df.subset.II, aes(x=HLA_Allele, y=RPKM)) + geom_boxplot(notch=FALSE, outlier.shape = NA, width=0.4) + geom_jitter(position=position_jitter(0.15),aes(color=`Tumor Type`), size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=14)) + ylab("RPKM") + scale_color_manual(values = my_pal) + xlab("Class II HLA Alleles")
plot(my_plot)

# same but x axis is sample
ggplot(df.tidy, aes(x=Samples, y=RPKM, color=Tumor)) + geom_point(aes(shape=`HLA Class`), size=3) + facet_wrap(~ Batch, scales = "free") + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=16)) + scale_color_manual(values = my_pal)

# dot plot with x axis is tumor type and legend is HLA allele
ggplot(df.tidy, aes(x=`Tumor Type`, y=RPKM)) + geom_point(position=position_jitter(0.15),aes(color=Tumor), size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=14)) + ylab("HLA Alleles (RPKM)") + scale_color_manual(values = my_pal)
