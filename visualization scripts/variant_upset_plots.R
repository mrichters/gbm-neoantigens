# variant upset plots #
# input is 3 excel files with annotated variants #

library(readxl)
library(tidyverse)
library(ggsci)
library(plyr); library(dplyr)
library(GGally)
library(UpSetR)
library(gridExtra)
library(gtable)
library(scales)
library(foreach)
library(mgsub)
library(scales)

# functions
filter_variants <- function(variant_df) {
  filtered_df <- data.frame()
  for (row in 1:nrow(variant_df)) {
    caller_set <- unlist(strsplit(variant_df[row, "set"], "-"))
    if (length(caller_set) > 1) {
      filtered_df <- rbind(filtered_df, variant_df[row, ])
    }
  }
  return(filtered_df)
}

divide_by <- function(x,y) { return(x/y) }

add_tumor_type <- function(df_column) {
  type_column <- c()
  for ( item in df_column ) {
    if ( grepl("BrMET", item) ) {
      tumor_type <- "BrMET"
    } else if ( grepl("Re", item) ) {
      tumor_type <- "Recurrent GBM"
    } else {
      tumor_type <- "GBM"
    }
  }
  type_column <- c(type_column, tumor_type)
  return(type_column)
}

setwd("~/Desktop/GBM_plots")
# list 3 gbm batches for analysis
gbm_batches <- c("Batch1", "Batch2", "Batch3")
tumor_count = 0
mut_df <- c()
total_df <- c()
clonal_df <- c()
for ( batch in gbm_batches[3] ) {
  variant_path <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_annotated_variants.xlsx", sep="")
  samples <- excel_sheets(path = variant_path)
  tumors <- unique(sapply(samples, function(x) return(unlist(strsplit(x, '-'))[1])))
  # if ( batch == "Batch3") {
  #   tumors <- unique(sapply(mgsub(samples, c("19-", "_"), c("", "-")), function(v) return(unlist(strsplit(v, '-'))[1])))
  # } else {
  #   tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[2])))
  # }
  for ( tumor in tumors[1] ) {
    tumor_count = tumor_count + 1
    tumor_df <- c()
    sample_names <- c()
    tumor_samples <- samples[grep(tumor, samples)]
    for ( name in samples[grep(tumor, samples)] ) {
      # read in correct batch variant excel file
      df <- as.data.frame(read_excel(variant_path, sheet = name))
      # if want to filter results to only include those called by 2 variant callers
      #df <- filter_variants(df)
      # create standardized shorter sample names
      # if ( grepl("Re", name) == TRUE ) {
      #     if ( grepl("Primary", name) == TRUE ) {
      #         new_name <- "Re.GBM065-primary"
      #         tumor <- "GBM065"
      #     } else {
      #     new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
      #     tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
      #     }
      # } else if ( batch == "Batch3") { 
      #     new_name <- mgsub(name, c("19-", "Tumor", "_"), c("", "", "-"))
      #     tumor <- unlist(strsplit(new_name, '-'))[1]
      # } else if ( grepl("Re", name) == FALSE ) {
      #     new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
      #     tumor <- unlist(str_split(name, '-'))[2]
      # }
      new_name <- name
      sample_names <- c(sample_names, new_name)
      #assign_name <- gsub("-", "\\.", new_name)
      # to write final dfs to excel file
      #assign(assign_name, df)
      # subset df
      df.subset <- df %>% 
        select(CHROM, POS, REF, ALT, SYMBOL, Consequence) 
      df.subset.cbind <- cbind(df.subset, Sample = new_name, Value = 1)
      df.final <- df.subset.cbind %>%
        unite(Variant, CHROM, POS, REF, ALT, sep=':')
      tumor_df <- rbind(tumor_df, df.final)
    }
    df.spread <- tumor_df %>% 
      spread(Sample, Value)
    df.spread[is.na(df.spread)] <- 0
    # to create df for variant clonality plots (clonal, subclonal shared, private)
    max_clonal_num <- ncol(df.spread[, -c(1:3)])
    df.spread$"Clonality_Count" <- rowSums(df.spread[, -c(1:3)])
    for ( line in 1:nrow(df.spread) ) {
      count <- df.spread[line, "Clonality_Count"]
      if ( max_clonal_num == 2 ) {
        if ( count == 1 ) {
        df.spread[line, "Coverage"] <- "Subclonal Private"
        } else if ( count == 2 ) {
          df.spread[line, "Coverage"] <- "Clonal"
        }
      } else {
        if ( count == 1 ) {
          df.spread[line, "Coverage"] <- "Subclonal Private"
        } else if ( count %in% 2:(max_clonal_num-1)) {
          df.spread[line, "Coverage"] <- "Subclonal Shared"
        } else if ( count == max_clonal_num ) {
          df.spread[line, "Coverage"] <- "Clonal"
        }
      }
    }
    
    final_spread <- data.frame(Tumor = tumor, as.data.frame(table(df.spread$Coverage)))
    colnames(final_spread) <- c("Tumor", "Variant Type", "Variant Count")
    clonal_df <- rbind(clonal_df, final_spread)
    # two waterfall specific dfs
    # Tumor Variant SYMBOL Consequence
    total_df <- rbind(total_df, data.frame(Tumor = tumor, select(df.spread, Variant, SYMBOL, Consequence)))
    # Tumor Total_Variants
    mut_df <- rbind(mut_df, data.frame("Tumor" = tumor, "Total_Variants" = nrow(df.spread)))
    # for D3 plots
    set_symbols <- c("A", "B", "C", "D")
    names(sample_names) <- set_symbols[1:length(sample_names)]
    for (s in sample_names) {
      df.spread[s] <- ifelse(test = df.spread[s] == 1, names(sample_names)[which(sample_names==s)], NA)
    }
    df.spread1 <- df.spread %>% unite(set, unname(sample_names), sep = ",", na.rm=TRUE)
    df.spread1$set2 <- I(sapply(df.spread1$set, function(x) ifelse(test = grepl(",", x), yes = strsplit(x, ","), no = (list(x)))))
    df.spread1 <- df.spread1[,c("Variant", "set2")] 
    names(df.spread1) <- c("name", "set")
    df.spread1$r <- "8"
    # create upset plots
    #variant_upset <- upset(df.spread, order.by = "degree", mainbar.y.label = "Variant Intersections", sets.x.label = "Variants per Sample", sets = sort(sample_names, decreasing = TRUE), keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)
    #pdf(paste(tumor, ".variant.upset_unfiltered.pdf", sep = ""), width = 10, height = 8)
    #print(variant_upset)
    #dev.off()
  }

clonal_df$Tumor <- factor(clonal_df$Tumor, levels = unique(sort(as.character(clonal_df$Tumor))))
clonal_df <- clonal_df[order(clonal_df$Tumor), ]
clonal_df$`Variant Type` <- factor(clonal_df$`Variant Type`, levels = c("Subclonal Private", "Subclonal Shared", "Clonal"))
clonal_df <- clonal_df[order(clonal_df$`Variant Type`),]

# write filtered variant set to excel file
#WriteXLS::WriteXLS(c(Re.GBM065.1, Re.GBM065.2), "GBM065_filtered.xlsx", SheetNames = tumor_samples)

}

# create clonal
clonality <- ggplot(filter(clonal_df, Tumor %in% c("BrMET025", "BrMET028")), aes(x=Tumor, y=`Variant Count`, fill=`Variant Type`)) + geom_bar(stat="identity") + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=13)) + ylab("Variant Count") + xlab("Tumor") + scale_fill_manual(values = c("green", "blue", "red")) #+ scale_y_continuous(trans="log10", breaks=c(1e+2, 1e+4, 1e+6, 1e+8, 1e+10, 1e+12))#+ scale_y_continuous(trans="sqrt", breaks = c(100, 100, 1000, 5000)) 
pdf(paste("B3_clonality.pdf", sep = ""), width = 5, height = 7)
print(clonality)
dev.off()

clonal_df_spr <- spread(clonal_df, key=`Variant Type`, value=`Variant Count`)
clonal_df_spr[is.na(clonal_df_spr)] <- 0
WriteXLS::WriteXLS(clonal_df_spr, "variant_clonality.xlsx")  

# create bar plot with proportions of clonality instead of raw values
clonal_df_spr$Total <- rowSums(clonal_df_spr[, -1])
for ( column in colnames(clonal_df_spr[, c(2:4)]) ) {
  proportion_name <- paste(column, "Proportion", sep=" ")
  clonal_df_spr[proportion_name] <- mapply(divide_by, clonal_df_spr[column], clonal_df_spr["Total"])
}
clonal_df_spr <- gather(clonal_df_spr, key="Proportion Category", value="Proportion", colnames(clonal_df_spr[,6:8]))
clonal_df_spr$`Proportion Category` <- sapply(clonal_df_spr$`Proportion Category`, function(x) return(gsub(" Proportion", "",x)))
clonal_df_spr$`Proportion Category` <- factor(clonal_df_spr$`Proportion Category`, levels = c("Clonal", "Subclonal Shared", "Subclonal Private"))
clonal_df_spr <- clonal_df_spr[order(clonal_df_spr$`Proportion Category`),]

clonal_df_proportion <- spread(clonal_df_spr[,c("Tumor", "Proportion Category", "Proportion")], key=`Proportion Category`, value=Proportion)
WriteXLS::WriteXLS(clonal_df_proportion, "variant_clonality_proportion.xlsx")

clonal_df_spr$`Tumor Type` <- sapply(clonal_df_spr$Tumor, add_tumor_type)

# split violin plot
#ggplot(Tcell_norm_gather.ordered, aes(x=`Tumor Type`, y=Score, fill=Cell_Type)) + geom_violin(scale="width", trim=F) + scale_fill_manual(values = my_pal) + theme_bw() + theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=12)) + ylab("Normalized xCell score")

clonal_df_spr$`Tumor Type`[clonal_df_spr$`Tumor Type`=="Recurrent GBM"] <- "GBM"
# create clonal proportion bar plot
ggplot(clonal_df_spr, aes(x=Tumor, y=Proportion, fill=`Proportion Category`)) + geom_bar(stat="identity") + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=12), legend.position = "bottom") + scale_fill_manual(values = c("green", "blue", "red")) + scale_y_continuous(expand=c(0,0))

# create clonal proportion split boxplot
ggplot(clonal_df_spr, aes(x=`Tumor Type`, y=Proportion, fill=`Proportion Category`)) + geom_boxplot(position=position_dodge(0.75)) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.8) + theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=14), legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank()) + scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values = c("firebrick4",  "mediumslateblue", "forestgreen"))

# find shared variants within any tumor
shared_variants <- data.frame()
for ( row in 1:nrow(df.spread) ) {
  line <- as.vector(df.spread[row, c(4:9)])
  if ( line[1] == 0 && line[3] == 0 && sum(line) == 4 ) {
    add_line <- df.spread[row, c(1:3)] 
    shared_variants <- rbind(shared_variants, add_line)
  }
}

##### WATERFALL PLOTS FROM total_df / mut_df #####


# format df
levels(total_df$Tumor) <- unique(sort(as.character(total_df$Tumor)))
levels(mut_df$Tumor) <- unique(sort(as.character(mut_df$Tumor)))

# most frequently mutated genes
#total_df_genes <- filter(unique_total_df, variant_class != "synonymous_variant") #& variant_class != "protein_altering_variant")
#total_df
table_genes <- as.data.frame(table(total_df$SYMBOL))
table.sort <- table_genes[order(-table_genes$Freq), ]
top_genes_list <- as.vector(table.sort$Var1[1:50])

# input variants per tumor file
#mut_df.subset <- filter(mut_df, Tumor %in% total_df_genes$Tumor)
#mut_df.subset.sort <- mut_df.subset[order(mut_df.subset$sample), ]
#samples_df.subset <- filter(samples_df, sample %in% total_df_genes$sample)
#samples_df.subset.sort <- samples_df.subset[order(samples_df.subset$sample), ]
sample_ordering_waterfall <- unique(sort(as.character(mut_df$Tumor)))

waterfall(total_df, fileType = "Custom", variant_class_order = mutation_priority, mainXlabel = TRUE, mutBurden = mut_df.subset.sort, mainPalette=mutationColours, clinData = samples_df.subset.sort, clinLegCol=1, clinVarCol=c('Brain Metastasis'='#c2ed67', 'Primary GBM'='#E63A27', 'Recurrent GBM'='#e69127'), clinVarOrder=c('Brain Metastasis', 'Primary GBM', 'Recurrent GBM'), section_heights=c(1, 5, 0.7)) #, sampOrder = sample_ordering_waterfall)







