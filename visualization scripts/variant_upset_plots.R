# variant upset plots & clonality plots & waterfall plots #
# input is excel file with annotated variants #

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
library(GenVisR)
library(ggsignif)

# functions
divide_by <- function(x,y) { return(x / y) }

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

# set wd
setwd("~/Desktop/GBM\ Figures")

# variant samples / tumors
variant_path <- "~/Desktop/GBM/gbm_google_drive/merged_annotated_variants.xlsx"
samples <- excel_sheets(path = variant_path)
tumors <- unique(sapply(samples, function(x) return(unlist(strsplit(x, '-'))[1])))

variant_count_df <- c()
total_df <- c()
clonal_df <- c()
for ( tumor in tumors ) {
  tumor_df <- c()
  sample_names <- samples[grep(tumor, samples)]
  for ( name in samples[grep(tumor, samples)] ) {
    # read in correct batch variant excel file
    df <- as.data.frame(read_excel(variant_path, sheet = name, col_names=T))
    # if want to filter results to only include those called by 2 variant callers (NA)
    #df <- filter_variants(df)
    # subset df
    df.subset <- df %>% 
      select(CHROM, POS, REF, ALT, SYMBOL, Consequence) 
    df.subset.cbind <- cbind(df.subset, Sample = name, Value = 1)
    df.final <- df.subset.cbind %>%
      unite(Variant, CHROM, POS, REF, ALT, sep=':')
    tumor_df <- rbind(tumor_df, df.final)
  }
  df.spread <- tumor_df %>% 
    spread(Sample, Value)
  df.spread[is.na(df.spread)] <- 0
  # create upset plots
  #variant_upset <- upset(df.spread, order.by = "degree", mainbar.y.label = "Variant Intersections", sets.x.label = "Variants per Sample", sets = sort(sample_names, decreasing = TRUE), keep.order = TRUE, text.scale = 1.75, point.size = 4, line.size = 1.5)
  #pdf(paste(tumor, ".variant.upset_ts.pdf", sep = ""), width = 10, height = 8)
  #print(variant_upset)
  #dev.off()
  #2 waterfall specific dfs
  #Tumor Variant SYMBOL Consequence
  total_df <- rbind(total_df, data.frame(Tumor = tumor, select(df.spread, Variant, SYMBOL, Consequence)))
  #Tumor Total_Variants
  variant_count_df <- rbind(variant_count_df, data.frame("Tumor" = tumor, "Total_Variants" = nrow(df.spread)))
  #df for variant clonality plots (clonal, subclonal shared, private)
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
  clonal_df$Tumor <- factor(clonal_df$Tumor, levels = get_gbm_order(tumors))
  clonal_df <- clonal_df[order(clonal_df$Tumor), ]
  clonal_df$`Variant Type` <- factor(clonal_df$`Variant Type`, levels = c("Subclonal Private", "Subclonal Shared", "Clonal"))
  clonal_df <- clonal_df[order(clonal_df$`Variant Type`),]
}

#variant_count_df$`Tumor Type` <- add_tumor_type(variant_count_df$Tumor)
# create clonal bar plot (prefer python for this one)
clonality <- ggplot(clonal_df, aes(x=Tumor, y=`Variant Count`, fill=`Variant Type`)) + geom_bar(stat="identity") + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8), legend.position = "none") + ylab("Variant Count") + xlab("Tumor") + scale_fill_manual(values = c("forestgreen", "royalblue3", "firebrick3")) + scale_y_continuous(trans="sqrt", breaks = c(100, 100, 1000, 5000)) + coord_cartesian(expand = F)
pdf(paste("ts_variant_clonality_barplot.pdf", sep = ""), width = 6, height = 5)
print(clonality)
dev.off()

# create clonality df
clonal_df_spr <- spread(clonal_df, key=`Variant Type`, value=`Variant Count`)
clonal_df_spr[is.na(clonal_df_spr)] <- 0
WriteXLS::WriteXLS(clonal_df_spr, "variant_clonality.xlsx")  

# create bar plot with proportions of clonality instead of raw values
clonal_df_spr$Total <- rowSums(clonal_df_spr[, -1])
for ( column in colnames(clonal_df_spr[, c(2:4)]) ) {
  proportion_name <- paste(column, "Proportion", sep=" ")
  clonal_df_spr[proportion_name] <- mapply(divide_by, clonal_df_spr[column], clonal_df_spr["Total"])
}
# proportion 
clonal_df_spr <- gather(clonal_df_spr, key="Proportion Category", value="Proportion", colnames(clonal_df_spr[,6:8]))
clonal_df_spr$`Proportion Category` <- sapply(clonal_df_spr$`Proportion Category`, function(x) return(gsub(" Proportion", "",x)))
clonal_df_spr$`Proportion Category` <- factor(clonal_df_spr$`Proportion Category`, levels = c("Clonal", "Subclonal Shared", "Subclonal Private"))
# final dfs
clonal_df_spr <- clonal_df_spr[order(clonal_df_spr$`Proportion Category`),]
clonal_df_proportion <- spread(clonal_df_spr[,c("Tumor", "Proportion Category", "Proportion")], key=`Proportion Category`, value=Proportion)

WriteXLS::WriteXLS(clonal_df_proportion, "variant_clonality_proportion.xlsx")

clonal_df_spr$`Tumor Type` <- sapply(clonal_df_spr$Tumor, add_tumor_type)
clonal_df_proportion$`Tumor Type` <- sapply(clonal_df_proportion$Tumor, add_tumor_type)

# split violin plot
#ggplot(Tcell_norm_gather.ordered, aes(x=`Tumor Type`, y=Score, fill=Cell_Type)) + geom_violin(scale="width", trim=F) + scale_fill_manual(values = my_pal) + theme_bw() + theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_rect(color = "black", fill=NA, size=1), text = element_text(size=12)) + ylab("Normalized xCell score")

clonal_df_spr$`Tumor Type`[clonal_df_spr$`Tumor Type`=="Recurrent GBM"] <- "GBM"

sig_test <- function(df, variable) {
  x <- filter(clonal_df_proportion, `Tumor Type` == "BrMET")[variable]
  y <- filter(clonal_df_proportion, `Tumor Type` %in% c("GBM", "Recurrent GBM"))[variable]
  htest <- t.test(x,y)$p.value #using default parameters
  return(c(variable, htest))
}

for (i in c("Clonal", "Subclonal Shared", "Subclonal Private")) {
  print(sig_test(clonal_df_proportion, i))
}

# create clonal proportion bar plot
#ggplot(clonal_df_spr, aes(x=Tumor, y=Proportion, fill=`Proportion Category`)) + geom_bar(stat="identity") + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90")) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=12), legend.position = "bottom") + scale_fill_manual(values = c("green", "blue", "red")) + scale_y_continuous(expand=c(0,0))

# create clonal proportion split boxplot
clonality_boxplot <- ggplot(clonal_df_spr, aes(x=`Tumor Type`, y=Proportion, fill=`Proportion Category`)) + geom_boxplot(position=position_dodge(0.75)) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize=0.5) + theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black", size=1), panel.grid.major = element_line(colour = "grey90"), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=8), legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0,1.15), breaks = c(0,0.25,0.5,0.75,1)) + scale_fill_manual(values = c("firebrick3",  "royalblue3", "forestgreen")) + geom_signif(y_position=c(1, 1.05, 1.1), xmin=c(0.8,1, 1.2), xmax=c(1.7, 2, 2.2),annotation=c("***", "*", "*"), tip_length=0.01)

clonality + clonality_boxplot + plot_layout(widths = c(2,1.5))

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

# mut burden df
#levels(variant_count_df$Tumor) <- get_gbm_order(variant_count_df$sample)
names(variant_count_df) <- c("sample", "mut_burden")
# sample order
sample_order_wf <- get_gbm_order(variant_count_df$sample)
# custom df
levels(total_df$Tumor) <- unique(sort(as.character(total_df$Tumor)))
total_df_wf <- total_df[,-2]
names(total_df_wf) <- c("sample", "gene", "variant_class")
total_df_wf$variant_class <- sapply(total_df_wf$variant_class, function(x) unlist(strsplit(x, "&"))[1])
total_df_wf$variant_class <- sapply(total_df_wf$variant_class, function(x) gsub("_", " ", x))
total_df_wf$variant_class <- sapply(total_df_wf$variant_class, function(x) str_to_title(x))
# collapse to only count 1 mutation per gene per tumor
#collapsed_df <- total_df_wf[,c(1,2)] %>% filter(sample != "GBM065.Re") %>% arrange(sample, gene) %>% distinct()
#gene_table <- as.data.frame(table(collapsed_df$gene))
#names(gene_table) <- c("gene", "frequency")

# get_genes <- function(df, tumor_type) {
#   df <- filter(df, grepl(tumor_type, sample))
#   df <- filter(df, variant_class %in% c("missense_variant", "frameshift_variant", "inframe_deletion", "start_lost", "stop_gained", "stop_lost", "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"))
#   collapsed_df <- df[,c(1,2)] %>% arrange(sample, gene) %>% distinct()
#   gene_table <- as.data.frame(table(collapsed_df$gene))
#   names(gene_table) <- c("gene", "frequency")
#   top_genes <- gene_table %>% filter(frequency >= 3) %>% arrange(desc(frequency))
#   gene_list <- top_genes$gene
#   return(top_genes)
# }
# #%>% filter(sample != "GBM065.Re")
# brmet <- get_genes(total_df_wf, "BrMET") 
# gbm_w_hypermutator <- get_genes(total_df_wf, "GBM") 
clin_data <- read_excel("~/Desktop/GBM Figures/clinical_data.xlsx")

 format_waterfall <- function(df, tumor_type, clin_col, clin_order) {  
  x <- filter(df, grepl(tumor_type, sample) & variant_class %in% c("Missense Variant", "Frameshift Variant", "Inframe Deletion", "Start Lost", "Stop Gained", "Stop Lost", "Splice Region Variant", "Splice Donor Variant", "Splice Acceptor Variant"))
  sample_order <- get_gbm_order(x$sample)
  x$sample <- factor(x$sample, levels = sample_order)
  #x$variant_class <- sapply(x$variant_class, function(x) unlist(strsplit(x, "&"))[1])
  #x <- x[order(x$sample),]
  variant_order <- unique(x$variant_class)
  mut_df <- filter(variant_count_df, grepl(tumor_type, sample))
  #mut_df$sample <- factor(mut_df$sample, levels = sample_order)
  #mut_df <- mut_df[order(mut_df$sample),]
  plot_waterfall(x, mut_df, sample_order, variant_order, clin_col, clin_order)
  #return(wf)
}

plot_waterfall <- function(df, mut_df, sample_order, variant_order, clin_col, clin_order) {
  GenVisR::waterfall(df, fileType = "Custom", variant_class_order = variant_order, mainXlabel = TRUE, mutBurden = mut_df, mainDropMut = TRUE, sampOrder = sample_order, rmvSilent = TRUE, mutBurdenLayer = list(scale_y_continuous(trans="log10")), maxGenes = 25, mainPalette = c("Stop Gained"='#4f00A8', "Frameshift Variant"='#A80100', "Inframe Deletion"='#ff9b34', "Missense Variant"='#009933', "Splice Region Variant" = "palevioletred"), clinDat = clin_data, clinVarCol = clin_col, clinLegCol = 2, clinVarOrder = clin_order)
}
# plotGenes = gene_list
gbm_waterf <- format_waterfall(total_df_wf, "GBM", c("Yes" = "lightgreen", "No" = "khaki1", "GBM" = "#3C5488FF", "Recurrent GBM" = "#F39B7FFF"), c("Yes", "No", "GBM", "Recurrent GBM")) #, gbm_w_hypermutator$gene)
brmet_wf <- format_waterfall(total_df_wf, "BrMET", c("Yes" = "lightgreen", "No" = "khaki1", "Breast Cancer" = "#E64B35FF", "Melanoma" = "#7E6148FF", "NSCLC" = "#00A087FF", "SCLC" = "#4DBBD5FF"), c("Yes", "No", "Breast Cancer", "Melanoma", "NSCLC", "SCLC")) #, brmet$gene)


#variant_class_order = c("missense_variant", "frameshift_variant", "inframe_deletion", "start_lost", "stop_gained", "stop_lost", "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"), mainPalette = c("stop_gained"='#4f00A8', "frameshift_variant"='#A80100', "inframe_deletion"='#ff9b34', "missense_variant"='#009933', "splice_region_variant" = "palevioletred")

#GenVisR::waterfall(gbm_df, fileType = "Custom", variant_class_order = c("missense_variant", "frameshift_variant", "inframe_deletion", "start_lost", "stop_gained", "stop_lost", "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant"), mainPalette = c("stop_gained"='#4f00A8', "frameshift_variant"='#A80100', "inframe_deletion"='#ff9b34', "missense_variant"='#009933', "splice_region_variant" = "palevioletred"), mainRecurCutoff = 0.05, maxGenes = 30, mainXlabel = TRUE, rmvSilent = TRUE, mutBurden = gbm_mut, mainDropMut = TRUE, sampOrder = sample_ordering_waterfall)#, clinLegCol=1, clinVarCol=c('Brain Metastasis'='#c2ed67', 'Primary GBM'='#E63A27', 'Recurrent GBM'='#e69127'), clinVarOrder=c('Brain Metastasis', 'Primary GBM', 'Recurrent GBM'))

#GenVisR::waterfall(brmet, fileType = "Custom", variant_class_order = mutation_priority, mainXlabel = TRUE, mutBurden = mut_df.subset.sort, mainPalette=mutationColours, clinData = samples_df.subset.sort, clinLegCol=1, clinVarCol=c('Brain Metastasis'='#c2ed67', 'Primary GBM'='#E63A27', 'Recurrent GBM'='#e69127'), clinVarOrder=c('Brain Metastasis', 'Primary GBM', 'Recurrent GBM'), section_heights=c(1, 5, 0.7)) #, sampOrder = sample_ordering_waterfall)


# for D3 plots
# set_symbols <- c("A", "B", "C", "D")
# names(sample_names) <- set_symbols[1:length(sample_names)]
# for (s in sample_names) {
#   df.spread[s] <- ifelse(test = df.spread[s] == 1, names(sample_names)[which(sample_names==s)], NA)
# }
# df.spread1 <- df.spread %>% unite(set, unname(sample_names), sep = ",", na.rm=TRUE)
# df.spread1$set2 <- I(sapply(df.spread1$set, function(x) ifelse(test = grepl(",", x), yes = strsplit(x, ","), no = (list(x)))))
# df.spread1 <- df.spread1[,c("Variant", "set2")] 
# names(df.spread1) <- c("name", "set")
# df.spread1$r <- "8"



