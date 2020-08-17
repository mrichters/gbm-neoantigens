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


# functions
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

# set wd
setwd("~/Desktop/GBM_plots")

# neoantigen upset plots #
# run for clinical / clinical no rna filter, classI / classII / total #

# list 3 gbm batches for analysis
gbm_batches <- c("Batch1", "Batch2", "Batch3")
clinical <- "clinical_no_rna_filter"

mut_df <- c()
total_df <- c()
clonal_df <- c()
for ( batch in gbm_batches ) {
  naI_path <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_pvacseq/", clinical, "/", batch, "_", clinical, "_classI_pvacseq.xlsx", sep="")
  naII_path <- paste("~/Desktop/GBM/gbm_google_drive/", batch, "_pvacseq/", clinical, "/", batch, "_", clinical, "_classII_pvacseq.xlsx", sep="")
  samples <- excel_sheets(path = naI_path)
  if ( batch == "Batch3") {
    tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(mgsub(v, c("19-", "_"), c("", "-")), "-"))[1])))
  } else {
    tumors <- unique(sapply(samples, function(v) return(unlist(strsplit(v, '-'))[2])))
  }
  for ( tumor in tumors ) {
    tumor_df <- c()
    sample_names <- c()
    for ( name in samples[grep(tumor, samples)] ) {
      # combine both for total neoantigens
      #df <- as.data.frame(rbind(read_excel(naI_path, sheet = name, cell_cols("A:AI")), read_excel(naII_path, sheet = name, cell_cols("A:AI"))))
      # read in correct batch neoantigen (naI or naII) excel file
      df <- as.data.frame(read_excel(naI_path, sheet = name, cell_cols("A:AI")))
      if ( nrow(df) != 0 ) {
        # create standardized shorter sample names
        if ( grepl("Re", name) == TRUE ) {
          if ( grepl("primary", name) == TRUE ) {
            new_name <- "Re.GBM065-primary"
            tumor <- "GBM065"
          } else {
            new_name <-  paste(paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = "."), str_sub(name, -1, -1), sep = '-')
            tumor <- paste("Re", str_sub(unlist(str_split(name, '-'))[2], 1, 6), sep = ".")
          }
        } else if ( batch == "Batch3") { 
          new_name <- paste(tumor, str_sub(name, -1, -1), sep = '-')
        } else if ( grepl("Re", name) == FALSE ) {
          new_name <- paste(unlist(str_split(name, '-'))[2], str_sub(name, -1, -1), sep = '-')
          tumor <- unlist(str_split(name, '-'))[2]
        }
        sample_names <- c(sample_names, new_name)
        # subset data frame to find unique neoantigens
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
    }
  }
}
clonal_df$Tumor <- factor(clonal_df$Tumor, levels = unique(sort(as.character(clonal_df$Tumor))))
clonal_df <- clonal_df[order(clonal_df$Tumor), ]
clonal_df$`Variant Type` <- factor(clonal_df$`Variant Type`, levels = c("Subclonal Private", "Subclonal Shared", "Clonal"))
clonal_df <- clonal_df[order(clonal_df$`Variant Type`),]
    # create upset plots
    #na_upset <- upset(df.spread, order.by = "degree", mainbar.y.label = "Neoantigen Intersections", sets.x.label = "ClassI Neoantigens per Sample", sets = sort(sample_names, decreasing = TRUE), keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)
    #pdf(paste(tumor, ".classI_neoantigen.clinical.nornafilter.upset.pdf", sep = ""), width = 12, height = 8)
    #print(na_upset)
    #dev.off()

# create clonal
ggplot(clonal_df, aes(x=Tumor, y=`Variant Count`, fill=`Variant Type`)) + geom_bar(stat="identity") + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=13)) + ylab("log10(Variant Count)") + xlab("Tumor") + scale_fill_manual(values = c("green", "blue", "red")) #+ scale_y_continuous(trans="log10", breaks=c(1e+2, 1e+4, 1e+6, 1e+8, 1e+10, 1e+12))#+ scale_y_continuous(trans="sqrt", breaks = c(100, 100, 1000, 5000)) 

clonal_df_spr <- spread(clonal_df, key=`Variant Type`, value=`Variant Count`)
clonal_df_spr[is.na(clonal_df_spr)] <- 0
WriteXLS::WriteXLS(clonal_df_spr, "neoantigen_clonality.xlsx")

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
WriteXLS::WriteXLS(clonal_df_proportion, "neoantigen_clonality_proportion.xlsx")

clonal_df_spr$`Tumor Type` <- sapply(clonal_df_spr$Tumor, add_tumor_type)
clonal_df_proportion$`Tumor Type` <- sapply(clonal_df_proportion$Tumor, add_tumor_type)

x <- filter(clonal_df_proportion, `Tumor Type` == "BrMET")$`Clonal`
y <- filter(clonal_df_proportion, `Tumor Type` %in% c("GBM", "Recurrent GBM"))$`Clonal`
t.test(x,y,var.equal=T)

clonal_df_spr$`Tumor Type`[clonal_df_spr$`Tumor Type`=="Recurrent GBM"] <- "GBM"
# create clonal proportion split boxplot
ggplot(clonal_df_spr, aes(x=`Tumor Type`, y=Proportion, fill=`Proportion Category`)) + geom_boxplot(position=position_dodge(0.75)) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.8) + theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=14), legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank()) + scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values = c("firebrick4",  "mediumslateblue", "forestgreen"))
