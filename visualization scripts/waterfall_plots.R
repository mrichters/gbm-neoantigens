# 8.29.19

# import packages
library(GenVisR)
library(readxl)
library(tidyverse)
library(ggsci)
library(plyr); library(dplyr)
library(GGally)
library(UpSetR)
library(gridExtra)
library(gtable)
library(scales)
library(mgsub)

#example_data <- read.table("~/Downloads/BKM120_Mutation_Data.tsv", header = TRUE, sep = '\t')

# Reformat the mutation data for waterfall()
# for y axis range!!
#+ scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))

#premade gene list
#genes_list <- c("TERT", "EGFR", "TP53", "PIK3CA", "PTEN", "MTOR")

#finding top 25 mutated genes in all my samples

#set wd
setwd("~/Desktop/GBM_plots")

gbm_batches <- c("Batch1", "Batch2", "Batch3")

# total_df == mutated genes per sample
# unique_total_df == mutated genes per tumor
total_df <- c()
#unique_total_df <- c()
mut_df <- c()
samples_df <- c()
for ( batch in gbm_batches ) {
  variant_file <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_annotated_variants.xlsx", sep = "")
  samples <- excel_sheets(path=variant_file)
  print(samples)
  if ( batch == "Batch3") {
    tumors <- unique(sapply(mgsub(samples, "19-", ""), function(v) return(unlist(strsplit(v, '_'))[1])))
  } else {
    tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[2])))
  }
  print(tumors)
  for ( tumor in tumors ) {
    unique_tumor_df <- c()
    total_tumor_df <- c()
    for ( name in samples[grep(tumor, samples)] ) {
      print(name)
      gbm_data <- as.data.frame(read_excel(variant_file, sheet = name))
      # creating shorter names for the plot
      if ( grepl("Re", name) == FALSE ) {
        if ( grepl("^19-", name) ) {
          new_name <- mgsub(name, c("19-", "_", "Tumor"), c("", "-", ""))
          #tumor <- unlist(str_split(name, '-'))[3]
        } else {
          new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
          #tumor <- unlist(str_split(name, '-'))[2]
        }
        if ( grepl("BrMET", name) ) {
          tumor_type <- "Brain Metastasis"
        } else {
          tumor_type <- "Primary GBM"
        }
      } else { 
        new_name <- paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
        #tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
        tumor_type <- "Recurrent GBM"
      }
      
      # format final df
      gbm_mutation_data <- gbm_data[,c("SYMBOL", "Consequence", "HGVSc")]
      gbm_mutation_data$Consequence <- unlist(lapply(strsplit(gbm_mutation_data$Consequence, "&"), `[[`, 1))
      gbm_mutation_data[is.na(gbm_mutation_data)] <- 0
      for ( row in 1:nrow(gbm_mutation_data) ) {
        if ( gbm_mutation_data$HGVSc[row] != 0 ) {
          gbm_mutation_data$HGVSc[row] <- unlist(lapply(strsplit(gbm_mutation_data$HGVSc[row], ":"), tail, 1))
        }
      }
      gbm_mutation_data <- cbind("sample" = new_name, gbm_mutation_data)
      colnames(gbm_mutation_data) <- c("sample", "gene", "variant_class", "variant")
      gbm_mutation_data_filter <- gbm_mutation_data #filter(gbm_mutation_data, variant != 0)
      # add data to dfs
      total_tumor_df <- rbind(total_tumor_df, gbm_mutation_data_filter)
      mut_df_line <- cbind(as.data.frame(new_name), as.data.frame(nrow(gbm_mutation_data_filter)))
      colnames(mut_df_line) <- c("sample", "mut_burden")
      mut_df <- rbind(mut_df, mut_df_line)
      samples_line <- data.frame(new_name, "Tumor Type", tumor_type)
      colnames(samples_line) <- c("sample", "variable", "value")
      samples_df <- rbind(samples_df, samples_line)
    }
    # create unique_tumor_df:
    #united_variants <- unite(total_tumor_df[, c(2,3,4)], "variants", sep=":")
    #uniqued_variants <- data.frame("unique.variants" = unique(sort(united_variants$variants)))
    #uniqued_variants.separate <- separate(uniqued_variants, unique.variants, c("gene", "variant_class", "amino.acid.change"), sep = ":")
    # add to unique_total_df
    #unique_total_df <- rbind(unique_total_df, data.frame("Tumor" = tumor, uniqued_variants.separate))
    # add to total_df
    total_df <- rbind(total_df, total_tumor_df)
  }
}
# total df / unique total df created!!

########## ORIGINAL BY SAMPLE CODE ###########
sample_order <- sort(as.vector(mut_df$sample))
# change sample order
total_df$sample <- as.character(total_df$sample)
total_df.sort <- total_df[order(total_df$sample),]
mut_df$sample <- as.character(mut_df$sample)
mut_df.sort <- mut_df[order(mut_df$sample), ]
# most frequently mutated genes
# new better way of picking genes
#total_df_genes <- filter(total_df, variant_class != "synonymous_variant") #& variant_class != "protein_altering_variant")
#table_genes <- as.data.frame(table(total_df_genes$gene))
#table.sort <- table_genes[order(-table_genes$Freq), ]
#top_genes_list <- as.vector(table.sort$Var1[1:50])
consequences <- c("stop_gained", "frameshift_variant", "inframe_deletion", "missense_variant", "splice_region_variant", "protein_altering_variant")
# original way of picking genes
#table_genes_df <- as.data.frame(table(total_df$gene))
#table.sort_df <- table_genes[order(-table_genes$Freq), ]
#top_genes_list_df <- as.vector(table.sort$Var1[1:50])
total_df_genes <- filter(total_df, variant_class %in% consequences)
#mut_df.subset <- filter(mut_df, sample %in% total_df_genes$sample)
mut_df.sort <- mut_df[order(mut_df$sample), ]
#samples_df.subset <- filter(samples_df, sample %in% total_df_genes$sample)
samples_df.sort <- samples_df[order(samples_df$sample), ]
sample_ordering_waterfall <- sort(as.vector(mut_df$sample))

new_mut_names <- sapply(total_df_genes$variant_class, function(x) return(gsub('_', ' ', paste(toupper(substring(x, 1,1)), substring(x, 2), sep=""))))
total_df_genes$variant_class <- new_mut_names

#toupper(substring(c, 1,1))

# Create a vector to save mutation priority order for plotting
#mutation_priority <- as.character(unique(total_df_genes$variant_class))
mutation_priority <- c("Stop gained", "Frameshift variant", "Inframe deletion", "Missense variant", "Splice region variant", "Protein altering variant")
mutationColours <- c("Stop gained"='#4f00A8', "Frameshift variant"='#A80100', "Inframe deletion"='#ff9b34', "Missense variant"='#009933', "Splice region variant"='#ca66ae', "Protein altering variant"='#A80079')

# Create the plot - all samples
waterfall(total_df_genes, fileType = "Custom", variant_class_order = mutation_priority, mainXlabel = FALSE, mutBurden = mut_df.sort, mainPalette=mutationColours, clinData = samples_df.sort, clinLegCol=1, clinVarCol=c('Brain Metastasis'='#c2ed67', 'Primary GBM'='#E63A27', 'Recurrent GBM'='#e69127'), clinVarOrder=c('Brain Metastasis', 'Primary GBM', 'Recurrent GBM'), section_heights=c(1, 5, 0.7), mainDropMut = TRUE, maxGenes=30) #, sampOrder = sample_ordering_waterfall)

## WITH NEW WATERFALL ##
total_df_new <- total_df[, c("sample", "gene", "variant_class")]
colnames(total_df_new) <- c("sample", "gene", "mutation")
hierarchy <- data.frame("mutation" = c("stop_gained", "frameshift_variant", "inframe_deletion", "missense_variant", "splice_region_variant", "protein_altering_variant"), "color" = c('#4f00A8','#A80100','#ff9b34','#009933','#ca66ae','#A80079'))
Waterfall(total_df_new, mutation = c("stop_gained", "frameshift_variant", "inframe_deletion", "missense_variant", "splice_region_variant", "protein_altering_variant"), mutationHierarchy = hierarchy)

#### CREATE THE BY SAMPLE PLOT BRMETS/GBMs ONLY
total_df_genes_alt <- total_df_genes[grepl("GBM", total_df_genes$sample), ]
#total_df_genes_alt <- total_df_genes_alt[!grepl("Re.GBM065", total_df_genes_alt$sample), ]
mut_df.subset <- filter(mut_df, sample %in% total_df_genes_alt$sample)
mut_df.subset.sort <- mut_df.subset[order(mut_df.subset$sample), ]
samples_df.subset <- filter(samples_df, sample %in% total_df_genes_alt$sample)
samples_df.subset.sort <- samples_df.subset[order(samples_df.subset$sample), ]
sample_ordering_waterfall <- sort(as.vector(mut_df.subset$sample))

# Create a vector to save mutation priority order for plotting
as.character(unique(total_df_genes_alt$variant_class))
#mutation_priority <- c("stop_gained", "frameshift_variant", "missense_variant")
mutation_priority <- c("stop_gained", "frameshift_variant", "inframe_deletion", "missense_variant", "splice_region_variant")#, "protein_altering_variant")

mutationColours <- c("stop_gained"='#4f00A8', "frameshift_variant"='#A80100', "inframe_deletion"='#ff9b34', "missense_variant"='#009933', "splice_region_variant"='#ca66ae')#, "protein_altering_variant"='#A80079')

waterfall(total_df_genes_alt, fileType = "Custom", variant_class_order = mutation_priority, mainXlabel = TRUE, mutBurden = mut_df.subset.sort, mainPalette=mutationColours, clinData = samples_df.subset.sort, clinLegCol=1, clinVarCol=c('Brain Metastasis'='#c2ed67', 'Primary GBM'='#E63A27', 'Recurrent GBM'='#e69127'), clinVarOrder=c('Brain Metastasis', 'Primary GBM', 'Recurrent GBM'), section_heights=c(1, 5, 0.7)) #, sampOrder = sample_ordering_waterfall)

############ NEW BY TUMOR CODE #############

# format df
levels(unique_total_df$Tumor) <- unique(sort(as.character(unique_total_df$Tumor)))

# most frequently mutated genes
total_df_genes <- filter(unique_total_df, variant_class != "synonymous_variant") #& variant_class != "protein_altering_variant")
table_genes <- as.data.frame(table(total_df_genes$gene))
table.sort <- table_genes[order(-table_genes$Freq), ]
top_genes_list <- as.vector(table.sort$Var1[1:50])

# input variants per tumor file
mut_df <- as.data.frame(read_tsv("variant_count_per_tumor.tsv"))
mut_df <- tumor_variant_counts[order(tumor_variant_counts$Tumor), ]
mut_df.subset <- filter(mut_df, strsplit(sample, "-")[[1]][1] %in% total_df_genes$Tumor)
mut_df.subset.sort <- mut_df.subset[order(mut_df.subset$sample), ]
samples_df.subset <- filter(samples_df, sample %in% total_df_genes$sample)
samples_df.subset.sort <- samples_df.subset[order(samples_df.subset$sample), ]
sample_ordering_waterfall <- sort(as.vector(mut_df.subset$sample))

# most frequently mutated genes
table_genes <- as.data.frame(table(total_df$gene))
table.sort <- table_genes[order(-table_genes$Freq), ]
top_genes_list <- as.vector(table.sort$Var1[1:50])
total_df_genes <- filter(total_df, gene %in% top_genes_list & variant_class != "synonymous_variant" & variant_class != "protein_altering_variant")
####
mut_df.subset <- filter(mut_df, sample %in% total_df_genes$sample)
mut_df.subset.sort <- mut_df.subset[order(mut_df.subset$sample), ]
samples_df.subset <- filter(samples_df, sample %in% total_df_genes$sample)
samples_df.subset.sort <- samples_df.subset[order(samples_df.subset$sample), ]
sample_ordering_waterfall <- sort(as.vector(mut_df.subset$sample))