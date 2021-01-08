# 1.28.2020

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
library(foreach)
library(mgsub)
library(BSDA)
library(ggsignif)


# functions
divide_by <- function(x,y) { return(x / y) }

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
    type_column <- c(type_column, tumor_type)
  }
  return(type_column)
}

get_gbm_order <- function(vector) {
  re_samples <- unique(vector[grepl("Re", vector)])
  not_re_samples <- unique(vector[!grepl("Re", vector)])
  sample_order <- c(not_re_samples, re_samples)
  tumor_order <- unique(sapply(sample_order, function(x) return(substr(x, 1, nchar(x)-2))))
  order_list <- list(sample_order, tumor_order)
  return(sample_order)
}

# neoag samples / tumors
neoag_path <- "~/Desktop/GBM/gbm_google_drive/merged_neoantigen_classII_filtered.xlsx"
samples <- excel_sheets(path = neoag_path)
tumors <- unique(sapply(samples, function(x) return(unlist(strsplit(x, '-'))[1])))

# set wd
setwd("~/Desktop/GBM Figures/")

mut_df <- c()
total_df <- c()
clonal_df <- c()
for ( tumor in tumors ) {
  tumor_df <- c()
  sample_names <- c()
  for ( name in samples[grep(tumor, samples)] ) {
    df <- as.data.frame(read_excel(neoag_path, sheet = name))
    new_name <- name
    if ( nrow(df) != 0 ) {
      # subset data frame to find unique neoantigens
      sample_names <- c(sample_names, name)
      df.subset <- df %>% 
        dplyr::select(Mutation, `Protein Position`, `Gene Name`, `MT Epitope Seq`)
      df.subset.cbind <- cbind(df.subset, Sample = new_name, Value = 1)
      df.final <- df.subset.cbind %>%
        unite(Neoantigen, Mutation, `Protein Position`, `Gene Name`, `MT Epitope Seq`, sep='_')
      tumor_df <- rbind(tumor_df, df.final)
    }
  }
  df.spread <- unique(tumor_df) %>% 
    spread(Sample, Value)
  df.spread[is.na(df.spread)] <- 0
  #if (ncol(df.spread) > 2 ) {
    #GBM065_na_upset <- upset(df.spread, order.by = "degree", mainbar.y.label = "Neoantigen Intersections", sets.x.label = "ClassII Neoantigens per Sample", sets = sort(sample_names, decreasing = TRUE), keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)
    #pdf(paste(tumor, ".classII_neoantigen.upset.pdf", sep = ""), width = 13, height = 8)
    #print(na_upset)
    #dev.off()
  # Tumor total neoantigens
  mut_df <- rbind(mut_df, data.frame("Tumor" = tumor, "Neoantigens" = nrow(df.spread)))
  # to create df for variant clonality plots (clonal, subclonal shared, private)
  max_clonal_num <- ncol(df.spread[, -1])
  if (is.null(max_clonal_num) == F) {
    df.spread$"Clonality_Count" <- rowSums(df.spread[, -1])
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
    clonal_df$Tumor <- factor(clonal_df$Tumor, levels = unique(sort(as.character(clonal_df$Tumor))))
    clonal_df <- clonal_df[order(clonal_df$Tumor), ]
    clonal_df$`Variant Type` <- factor(clonal_df$`Variant Type`, levels = c("Subclonal Private", "Subclonal Shared", "Clonal"))
    clonal_df <- clonal_df[order(clonal_df$`Variant Type`),]
  }
}

clonal_df$Tumor <- factor(clonal_df$Tumor, levels = get_gbm_order(tumors))
# create clonal
na_clonality <- ggplot(clonal_df, aes(x=Tumor, y=`Variant Count`, fill=`Variant Type`)) + geom_bar(stat="identity") + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8), axis.title.x = element_blank(), legend.position = "none") + ylab("Neoantigen Count") + xlab("Tumor") + scale_fill_manual(values = c("forestgreen", "royalblue3", "firebrick3")) + scale_y_continuous(trans="sqrt", breaks = c(100,500,1000,2000,4000,6000)) + coord_cartesian(expand = F)

clonal_df_spr <- spread(clonal_df, key=`Variant Type`, value=`Variant Count`)
clonal_df_spr[is.na(clonal_df_spr)] <- 0
WriteXLS::WriteXLS(clonal_df_spr, "neoantigen_clonality_bar_II.xlsx")

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
WriteXLS::WriteXLS(clonal_df_proportion, "neoantigen_clonality_bar_II.xlsx")

clonal_df_spr$`Tumor Type` <- sapply(clonal_df_spr$Tumor, add_tumor_type)
clonal_df_proportion$`Tumor Type` <- sapply(clonal_df_proportion$Tumor, add_tumor_type)

sig_test <- function(df, variable) {
  x <- filter(clonal_df_proportion, `Tumor Type` == "BrMET")[variable]
  y <- filter(clonal_df_proportion, `Tumor Type` %in% c("GBM", "Recurrent GBM"))[variable]
  htest <- t.test(x,y,var.equal=T)$p.value
  return(c(variable, htest))
}

for (i in c("Clonal", "Subclonal Shared", "Subclonal Private")) {
  print(sig_test(clonal_df_proportion, i))
}

clonal_df_spr$`Tumor Type`[clonal_df_spr$`Tumor Type`=="Recurrent GBM"] <- "GBM"
# create clonal proportion split boxplot
na_clonal_prop <- ggplot(clonal_df_spr, aes(x=`Tumor Type`, y=Proportion, fill=`Proportion Category`)) + geom_boxplot(position=position_dodge(0.75)) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.5) + theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8), axis.title.x = element_blank(), legend.position = "bottom") + scale_y_continuous(expand=c(0,0.01), limits = c(0, 1.05)) + scale_fill_manual(values = c("firebrick3",  "royalblue3", "forestgreen")) + geom_signif(y_position=c(0.9, 0.95, 1.02), xmin=c(0.8,1, 1.2), xmax=c(1.7, 2, 2.2),annotation=c("***", "N.S.", "**"), tip_length=0.01)
