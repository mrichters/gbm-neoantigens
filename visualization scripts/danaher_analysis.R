# 11.6.19
# perform danaher method for estimating immune infiltration and make a heatmap of the results

# import packages
library(tidyverse)
library(readxl)
library(openxlsx)
library(ComplexHeatmap)
library(ggsci)
library(rapportools)
library(wesanderson)
library(ggdendro)
library(gridExtra)
library(colorRamps)

# set wd
setwd("~/Desktop/GBM/immune_infiltration/")
## for expression: log2(fpkm+0.001) / log2(tpm+0.001)

# 9.21.20
######### KALLISTO GENE EXPRESSION EXCEL FILE TO MATRIX FOR XCELL AND CIBERSORT ###########
excel_file <- read_excel("kallisto_abundance_matrix_strandfix.xlsx")
samples <- names(excel_file[,2:ncol(excel_file)])
# filter kallisto file by danaher gene set
dgene <- read_excel("danaher_analysis/danaher_gene_set.xlsx")
excel_file <- filter(excel_file, gene_name %in% dgene$Gene)
# merge two files to include Cell Types
final_df <- merge(dgene, excel_file, by.x = "Gene", by.y = "gene_name")
# which gene not present in kallisto output?
not_included_sample <- filter(dgene, !Gene %in% excel_file$gene_name)
# log2 + 1 transform final_df
final_df_norm <- data.frame(final_df[,c(1:2)], sapply(final_df[samples], function(x) log2(x + 1)))
# get sample order
gbm_order <- get_gbm_order(samples)
# gather samples to average expression per cell type
final_gather <- gather(final_df_norm, key="Samples", value="Norm. Expression", -Gene, -Cell.Type)
# average expression
avg_cell_types <- final_gather %>% group_by(Cell.Type, Samples) %>% summarize(`Avg. Norm. Expression` = median(`Norm. Expression`))
# spread and create matrix
final_spread <- spread(avg_cell_types, key=Samples, value=`Avg. Norm. Expression`)
# matrix for ComplexHeatmap
final_spread_matrix <- as.matrix(column_to_rownames(final_spread, var = "Cell.Type"))
# row order (cell similarity grouping)
row_order <- c("B-cells", "CD45", "CD8 T cells", "Cytotoxic cells", "Exhausted CD8", "T-cells", "Th1 cells", "Treg", "DC", "Macrophages", "Mast cells", "Neutrophils", "NK cells", "NK CD56dim cells")
# create Heatmap!
Heatmap(final_spread_matrix, name = "log2(TPM + 1)", row_order = row_order, width = unit(28, "cm"), height = unit(8, "cm"))
# column_order = gsub('-', '\\.', gbm_order)



# 11.11.19
#### POLE PATIENT ALL TUMORS ABUNDANCE MATRIX ####
pole_file <- read_delim("matrix.abundance.tsv", '\t')
new_colnames <- replace(colnames(pole_file), colnames(pole_file)=="X1", "Genes")
colnames(pole_file) <- new_colnames
gene_list <- read_excel("danaher_analysis/danaher_gene_set.xlsx")
gene_set <- gene_list$Gene
pole_file.subset <- data.frame()
for ( row in 1:nrow(gene_list) ) {
  gene <- gene_list$Gene[row]
  category <- gene_list$`Cell Type`[row]
  gene_info <- pole_file[grep(paste("^", gene, "$", sep = "" ), pole_file$Genes, perl=TRUE), ]
  if ( is.empty(gene_info) == FALSE ) {
    new_line <- cbind(as.data.frame(category), as.data.frame(gene_info))
    pole_file.subset <- rbind(pole_file.subset, new_line)
  } else {
    print(gene)
  }
}
pole_file.gather <- gather(pole_file.subset, key = Sample, value = Abundance, -category, -Genes)
pole_file.gather$Sample <- replace(pole_file.gather$Sample, pole_file.gather$Sample=="GBM07-GBM27_TumorRNA_SF4", "GBM07-Escape")

samples <- unique(pole_file.gather$Sample)
# read in all sheets from excel file
abundance_df <- data.frame()
for ( name in samples ) {
  df <- filter(pole_file.gather, Sample == name)
  categories <- unique(df$category)
  for ( item in categories ) {
    values <- filter(df, category == item)$Abundance
    final_list <- c()
    for ( value in values ) {
      abundance <- log2(value + 1)
      final_list <- c(final_list, abundance)
    }
    abundance_value <- mean(final_list)
    print(abundance_value)
    new_line <- cbind(Samples = name, Cell_Category = item, Normalized_Expression = abundance_value)
    abundance_df <- rbind(abundance_df, new_line)
  }
}
####


# functions - get correct sample order
get_gbm_order <- function(vector) {
    re_samples <- unique(vector[grepl("Re", vector)])
    not_re_samples <- unique(vector[!grepl("Re", vector)])
    sample_order <- c(not_re_samples, re_samples)
    tumor_order <- unique(sapply(sample_order, function(x) return(substr(x, 1, nchar(x)-2))))
    order_list <- list(sample_order, tumor_order)
    return(sample_order)
}

# final_df_norm 
# samples
# get sample order
gbm_order <- get_gbm_order(samples)
# gather samples to plot
final_gather <- gather(final_df_norm, key="Samples", value="Norm. Expression", -Gene, -Cell.Type)
final_norm_matrix <- as.matrix(column_to_rownames(final_df_norm[,-2], var = "Gene"))

#function for factor to numeric: as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
# convert factor to numeric
#abundance_df$Normalized_Expression <- as.numeric(levels(abundance_df$Normalized_Expression))[abundance_df$Normalized_Expression]
# get correct sample order and order on samples
# gbm_order <- get_gbm_order(abundance_df)
#abundance_order <- factor_and_order(abundance_df, abundance_df$Samples, gbm_order[[1]])
# abundance_df$Samples <- factor(abundance_df$Samples, levels = gbm_order[[1]])
# abundance_order <- abundance_df[order(abundance_df$Samples),]

abundance_order <- filter(abundance_order, !Cell_Category %in% c("Exhausted CD8", "T-cells", "Cytotoxic cells", "NK CD56dim cells"))

# order on cell types
cell_order <- rev(c("CD45", "CD4 T cells", "CD8 T cells", "Th1 cells", "Treg", "PD-1", "PD-L1", "DC", "B-cells", "Macrophages", "Mast cells", "Neutrophils", "NK cells"))
abundance_order$Cell_Category <- factor(abundance_order$Cell_Category, levels = cell_order)
# for COMPLEX HEATMAP
# spread sample names, expression values
abundance_df.format <- spread(abundance_df, Samples, Normalized_Expression)#key = sample, value = abundance))
# change cell categories to row names
abundance_df.format <- column_to_rownames(abundance_df.format, loc=1)
# convert to matrix
abundance_matrix <- as.matrix(abundance_df.format)
# specify row/column order            
row_orders <- rownames(abundance_df.format)
column_orders <- excel_samples
colnames(abundance_order) <- c("Samples", "Cell_Category", "Cell Type Score")
 
#### GGTILE ####
# create palette
pal = wes_palette("Zissou1", 100, type = "continuous")

# base tile plot
ggplot(abundance_order, aes(x=Samples,  y=Cell_Category, fill=`Cell Type Score`)) + geom_tile() + theme(axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text = element_text(size=11), legend.position = "bottom", legend.key.width=unit(1,"cm"), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio=1) + scale_fill_gradientn(colours = pal)


## COMPLEX HEATMAP ##
# spread sample names, expression values
#abundance_df.format <- spread(abundance_df, Samples, Normalized_Expression)#key = sample, value = abundance))
# change cell categories to row names
#abundance_df.format <- column_to_rownames(abundance_df.format, loc=1)
# convert to matrix
#abundance_matrix <- as.matrix(abundance_df.format)
# specify row/column order            
row_orders <- rownames(final_df_norm)
#column_orders <- excel_samples

# specify heatmap parameters
col_fun = colorRamp2(c(1, 0.2, 0), c("red", "white", "blue"))
col_fun(seq(-3, 3))

# create Heatmap!
Heatmap(final_norm_matrix, name = "Enrichment Score", column_order = gsub('-', '\\.', gbm_order)) #width = unit(28, "cm"), height = unit(15, "cm"))

#Heatmap(total_matrix, row_order = c("CD4","CD8A", "CD8B","CD3D", "CD3E", "CD3G", "PDCD1", "CD274"), name = "log(TPM)", width = unit(31, "cm"), height = unit(12, "cm"))
#Heatmap(total_matrix, name = "log(TPM)", width = unit(31, "cm"), height = unit(25, "cm"), column_order = samples)
Heatmap(abundance_matrix, name = "log2(TPM)", width = unit(28, "cm"), height = unit(8, "cm"), column_order = gbm_order[[1]])
